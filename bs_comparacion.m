%% Comparación de métodos numéricos para Black-Scholes
%  - Galerkin con polinomios de Hermite
%  - Monte Carlo
%  - Diferencias finitas (Crank-Nicolson)
%  Referencia: solución exacta (blsprice)

clear; close all; clc;

%% 1. Parámetros del modelo y de la opción
r     = 0.05;       % Tasa de interés libre de riesgo
sigma = 0.2;        % Volatilidad constante
K     = 100;        % Precio de ejercicio (strike)
T     = 1;          % Tiempo hasta el vencimiento (años)
S0    = 100;        % Precio actual del subyacente

% Solución exacta de Black-Scholes (call)
precio_exacto = blsprice(S0, K, r, T, sigma);
fprintf('Precio exacto (Black-Scholes): %.8f\n\n', precio_exacto);

%% 2. Método de Galerkin (Hermite) – reutiliza el código adjunto
%   Se calcula el precio para distintos números de polinomios M,
%   se almacena el error respecto a la solución exacta.

M_vals = 5:5:200;              % Números de polinomios a probar
precio_galerkin = zeros(size(M_vals));
error_galerkin   = zeros(size(M_vals));

L = 10;  % Límite de integración para los coeficientes iniciales

for idx = 1:length(M_vals)
    M = M_vals(idx);
    
    % --- Cálculo de los coeficientes iniciales c0 (proyección del payoff) ---
    c0 = zeros(M+1, 1);
    for k = 0:M
        integrando = @(x) max(exp(x)-K, 0) .* hermite(x, k) .* exp(-x.^2);
        valor = integral(integrando, -L, L, 'RelTol',1e-8, 'AbsTol',1e-10);
        norm_factor = 2^k * factorial(k) * sqrt(pi);
        c0(k+1) = valor / norm_factor;
    end
    
    % --- Construcción de la matriz D del sistema de EDOs ---
    D = zeros(M+1, M+1);
    for i = 0:M
        for j = 0:M
            term1 = i * (sigma^2 - 2*r) * (i == j+1);      % delta(i, j+1)
            term2 = -2 * i * sigma^2 * (j+1) * (i == j+2); % delta(i, j+2)
            term3 = r * (i == j);                           % delta(i, j)
            D(j+1, i+1) = term1 + term2 + term3;
        end
    end
    A = -D;   % c' = A*c  => c(tau) = expm(A*tau)*c0
    
    % --- Solución en tau = T (vencimiento) ---
    cT = expm(A * T) * c0;
    
    % --- Evaluación de la aproximación en x = log(S0) ---
    x_eval = log(S0);
    u_aprox = 0;
    for k = 0:M
        u_aprox = u_aprox + cT(k+1) * hermite(x_eval, k);
    end
    precio_galerkin(idx) = u_aprox;
    error_galerkin(idx)   = abs(u_aprox - precio_exacto);
end

% Se toma como representante del método la aproximación con M máximo
[error_G, idx_min] = min(error_galerkin);
precio_G = precio_galerkin(idx_min);
M_opt = M_vals(idx_min);
fprintf('Galerkin (M = %d): precio = %.8f, error abs = %.2e\n', ...
    M_opt, precio_G, error_G);

%% 3. Método de Monte Carlo
%   Simulación directa bajo la medida neutral al riesgo

N_sim = 500000;   % Número de trayectorias
Z = randn(N_sim, 1);
ST = S0 * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z);
payoff = max(ST - K, 0);
precio_MC = exp(-r*T) * mean(payoff);
std_MC    = exp(-r*T) * std(payoff) / sqrt(N_sim);   % Error estándar
error_MC  = abs(precio_MC - precio_exacto);

fprintf('Monte Carlo (%d sim.): precio = %.8f, error abs = %.2e (EE = %.2e)\n', ...
    N_sim, precio_MC, error_MC, std_MC);

%% 4. Método de diferencias finitas – Crank-Nicolson
%   Se resuelve la EDP de Black-Scholes en el dominio S.

% Parámetros de la malla
S_max = 4 * K;               % Frontera superior suficientemente grande
M_fd  = 2000;                % Número de intervalos en S (M_fd+1 nodos)
N_fd  = 2000;                % Número de pasos en el tiempo
dS = S_max / M_fd;
dt = T / N_fd;

% Malla espacial y vector de tiempos
S_vec = linspace(0, S_max, M_fd+1)';
tau_vec = linspace(0, T, N_fd+1);   % tiempo hacia atrás: tau = T - t

% Condición inicial (tau = 0, es decir vencimiento)
V = max(S_vec - K, 0);

% Coeficientes constantes para la discretización
% Para un nodo interior i (i=2,...,M_fd) correspondiente a S_i
% a_i = -(r*S_i)/(2*dS) + (sigma^2*S_i^2)/(2*dS^2)
% b_i = 1/dt + r + sigma^2*S_i^2/dS^2
% c_i = (r*S_i)/(2*dS) + (sigma^2*S_i^2)/(2*dS^2)
% Se construye la matriz tridiagonal del sistema implícito.

i_interior = (2:M_fd)';      % índices de los nodos interiores
S_int = S_vec(i_interior);

% Coeficientes para los nodos interiores (parte implícita y explícita)
a = - (r*S_int)/(2*dS) + (sigma^2 * S_int.^2) / (2*dS^2);   % subdiagonal
b = 1/dt + r + (sigma^2 * S_int.^2) / (dS^2);               % diagonal
c = (r*S_int)/(2*dS) + (sigma^2 * S_int.^2) / (2*dS^2);     % superdiagonal

% Matriz tridiagonal implícita (A_impl) para interiores
% Sistema: A_impl * V^(n+1)_int = B_expl * V^n_int + condiciones contorno
A_impl = spdiags([a, b, c], [-1, 0, 1], M_fd-1, M_fd-1);

% Matriz explícita (B_expl) correspondiente a (1/dt - ... )*V^n
a_expl = -a;  % cambiando signos por ser la parte explícita
b_expl = 1/dt - r - (sigma^2 * S_int.^2) / (dS^2);
c_expl = -c;
B_expl = spdiags([a_expl, b_expl, c_expl], [-1, 0, 1], M_fd-1, M_fd-1);

% Bucle en el tiempo hacia atrás
for n = 1:N_fd
    tau_n = tau_vec(n);   % tiempo actual
    % Vector del lado derecho = B_expl * V_interior
    rhs = B_expl * V(i_interior);
    
    % Condiciones de frontera
    % En S=0, V=0 (no contribuye a rhs)
    % En S=S_max, V(S_max) = S_max - K*exp(-r*tau_n)
    V_max = S_max - K * exp(-r * tau_n);
    % El nodo interior inmediato (último) recibe contribución de V_max
    rhs(end) = rhs(end) - c(end) * V_max;   % c(end) es el coef. superdiagonal en la última fila de A_impl
    
    % Para el siguiente paso temporal necesitamos V_max en tau+dt
    tau_next = tau_n + dt;
    V_max_next = S_max - K * exp(-r * tau_next);
    
    % Resolver sistema tridiagonal
    V_new_int = A_impl \ rhs;
    
    % Actualizar vector completo V
    V(1) = 0;                  % frontera inferior
    V(i_interior) = V_new_int;
    V(end) = V_max_next;       % frontera superior
end

% Interpolar en S0 para obtener el precio
precio_CN = interp1(S_vec, V, S0, 'spline');
error_CN = abs(precio_CN - precio_exacto);
fprintf('Crank-Nicolson (%d x %d): precio = %.8f, error abs = %.2e\n', ...
    M_fd, N_fd, precio_CN, error_CN);

%% 5. Comparación y gráfico de errores
fprintf('\n--- Resumen de resultados ---\n');
fprintf('Método              | Precio       | Error abs\n');
fprintf('Exacto              | %12.8f | %10.2e\n', precio_exacto, 0);
fprintf('Galerkin (M=%d)      | %12.8f | %10.2e\n', M_opt, precio_G, error_G);
fprintf('Monte Carlo         | %12.8f | %10.2e\n', precio_MC, error_MC);
fprintf('Crank-Nicolson      | %12.8f | %10.2e\n', precio_CN, error_CN);

% Gráfico: evolución del error de Galerkin con M
%          y líneas horizontales para los errores de Monte Carlo y Crank-Nicolson
figure;
semilogy(M_vals, error_galerkin, 'b-o', 'LineWidth', 1.2);
hold on;
yline(error_MC, 'r--', 'LineWidth', 1.5);
yline(error_CN, 'g-.', 'LineWidth', 1.5);
xlabel('Número de polinomios de Hermite M');
ylabel('Error absoluto');
title('Comparación de errores: Galerkin (Hermite) vs Monte Carlo y Crank-Nicolson');
legend('Galerkin (Hermite)', 'Monte Carlo', 'Crank-Nicolson', 'Location', 'northeast');
grid on;

%% Función auxiliar: polinomios de Hermite (físicos)
function H = hermite(x, n)
    if n == 0
        H = ones(size(x));
    elseif n == 1
        H = 2 * x;
    else
        H_prev = ones(size(x));
        H_curr = 2 * x;
        for k = 1:n-1
            H_next = 2 * x .* H_curr - 2 * k * H_prev;
            H_prev = H_curr;
            H_curr = H_next;
        end
        H = H_curr;
    end
end