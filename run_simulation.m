function [SIMULATION_DATA, PARAMETERS] = run_simulation(PARAMETERS)
     
    % Prepare simulation data
    SIMULATION_DATA = {};
    SIMULATION_DATA.STW = zeros(PARAMETERS.SIMULATION.TOTAL_STEPS, 5); % Time, state, control, disturbance, epsilon
    x_stw = PARAMETERS.SIMULATION.INITIAL_STATE;
    epsilon_stw = 0;
    SIMULATION_DATA.STW_SAT = zeros(PARAMETERS.SIMULATION.TOTAL_STEPS, 5); % Time, state, control, disturbance, epsilon
    x_stw_sat = PARAMETERS.SIMULATION.INITIAL_STATE;
    epsilon_stw_sat = 0;

    SIMULATION_DATA.SNSTA_1 = zeros(PARAMETERS.SIMULATION.TOTAL_STEPS, 6); % Time, state, control, disturbance, epsilon, omega_c
    x_asnsta_1 = PARAMETERS.SIMULATION.INITIAL_STATE;
    epsilon_asnsta_1 = 0;
    est_x_asnsta_1 = x_asnsta_1;
    control_asnsta_1 = 0;

    SIMULATION_DATA.SNSTA_2 = zeros(PARAMETERS.SIMULATION.TOTAL_STEPS, 6); % Time, state, control, disturbance, epsilon, omega_c
    x_asnsta_2 = PARAMETERS.SIMULATION.INITIAL_STATE;
    epsilon_asnsta_2 = 0;
    est_x_asnsta_2 = x_asnsta_2;
    control_asnsta_2 = 0;

    % Run simulation
    simulation_time = 0.0;
    for simulation_step = 1:PARAMETERS.SIMULATION.TOTAL_STEPS
    
        if PARAMETERS.SIMULATION.SCENARIO == 1
            disturbance = 0;
        elseif PARAMETERS.SIMULATION.SCENARIO == 2
            disturbance = 2 + 0.6 * sin(simulation_time) + 0.4 *sin(5*simulation_time);
        elseif PARAMETERS.SIMULATION.SCENARIO == 3
            disturbance = 1 + 1*(sin(2.45*simulation_time) + exp(-4.85*abs(cos(3.13*simulation_time - 0.52))));
        else
            disturbance = 2;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------------------------------------------------------------------------------------------------% 
        %------------------------------------------------------ STW ----------------------------------------------------------% 
        %---------------------------------------------------------------------------------------------------------------------% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PARAMETERS.SIMULATION.NOISE_MODULE_DB > 0
            noise_x_stw = awgn(x_stw , PARAMETERS.SIMULATION.NOISE_MODULE_DB, "measured");   
        else
            noise_x_stw = x_stw;
        end 
        dot_epsilon_stw = -PARAMETERS.CONTROL.K2*sign(noise_x_stw);
        epsilon_stw = epsilon_stw + dot_epsilon_stw*PARAMETERS.SIMULATION.SAMPLING_TIME;       
        control_stw = -PARAMETERS.CONTROL.K1*sqrt(abs(noise_x_stw))*sign(noise_x_stw) + epsilon_stw;
        control_stw = function_sat(control_stw, PARAMETERS.CONTROL.MAX);
        SIMULATION_DATA.STW(simulation_step,:) = [simulation_time; x_stw; control_stw; disturbance; epsilon_stw];
        dot_x_stw = control_stw + disturbance;
        x_stw = x_stw + dot_x_stw*PARAMETERS.SIMULATION.SAMPLING_TIME;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------------------------------------------------------------------------------------------------% 
        %----------------------------------------------- STW/SATURATED -------------------------------------------------------% 
        %---------------------------------------------------------------------------------------------------------------------% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PARAMETERS.SIMULATION.NOISE_MODULE_DB > 0
            noise_x_stw_sat = awgn(x_stw_sat , PARAMETERS.SIMULATION.NOISE_MODULE_DB, "measured");   
        else
            noise_x_stw_sat = x_stw_sat;
        end
        dot_epsilon_stw_sat = -PARAMETERS.CONTROL.K2*sign(noise_x_stw_sat)-PARAMETERS.CONTROL.K3*epsilon_stw_sat;
        epsilon_stw_sat = epsilon_stw_sat + dot_epsilon_stw_sat*PARAMETERS.SIMULATION.SAMPLING_TIME;
        control_stw_sat = -PARAMETERS.CONTROL.K1*function_sat(sqrt(abs(noise_x_stw_sat)), PARAMETERS.CONTROL.EPSILON)*sign(noise_x_stw_sat) + epsilon_stw_sat;
        control_stw_sat = function_sat(control_stw_sat, PARAMETERS.CONTROL.MAX);
        SIMULATION_DATA.STW_SAT(simulation_step,:) = [simulation_time; x_stw_sat; control_stw_sat; disturbance; epsilon_stw_sat];
        dot_x_stw_sat = control_stw_sat + disturbance;
        x_stw_sat = x_stw_sat + dot_x_stw_sat*PARAMETERS.SIMULATION.SAMPLING_TIME;
   
        %---------------------------------------------------------------------------------------------------------------------% 
        %------------------------------------------------ SNSTA PARAMETERS ---------------------------------------------------% 
        %---------------------------------------------------------------------------------------------------------------------% 
        lambda = PARAMETERS.CONTROL.LAMBDA;
        gamma = PARAMETERS.CONTROL.GAMMA;
        p2 = (5*lambda^2+4)/(2*lambda^3);
        beta = 0.5*((1/(2*PARAMETERS.SIMULATION.SAMPLING_TIME)) - lambda)/gamma;
        kappa = (beta*2)/(p2);
        % disp('------------------');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------------------------------------------------------------------------------------------------% 
        %----------------------------------------------------- SNSTA_1 --------------------------------------------------------% 
        %---------------------------------------------------------------------------------------------------------------------% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PARAMETERS.SIMULATION.NOISE_MODULE_DB > 0
            noise_x_asnsta_1 = awgn(x_asnsta_1 , PARAMETERS.SIMULATION.NOISE_MODULE_DB, "measured");   
        else
            noise_x_asnsta_1 = x_asnsta_1;
        end
        z_asnsta_1 = noise_x_asnsta_1 - est_x_asnsta_1;
        dot_est_x_asnsta_1 = lambda*z_asnsta_1 + beta*tanh(gamma*z_asnsta_1) + epsilon_asnsta_1 + control_asnsta_1;
        dot_epsilon_asnsta_1 = ((lambda^2)/4.0)*z_asnsta_1 + kappa*tanh(gamma*z_asnsta_1);
        epsilon_asnsta_1 = epsilon_asnsta_1 + dot_epsilon_asnsta_1*PARAMETERS.SIMULATION.SAMPLING_TIME;
        est_x_asnsta_1 = est_x_asnsta_1 + dot_est_x_asnsta_1*PARAMETERS.SIMULATION.SAMPLING_TIME;
        d_est_x_asnsta_1 = lambda*z_asnsta_1 + beta*tanh(gamma*z_asnsta_1) + epsilon_asnsta_1;
        U = PARAMETERS.CONTROL.MAX - abs(d_est_x_asnsta_1);
        control_asnsta_1 = -U*tanh(PARAMETERS.CONTROL.CHI_1*noise_x_asnsta_1) - d_est_x_asnsta_1;
        control_asnsta_1 = function_sat(control_asnsta_1, PARAMETERS.CONTROL.MAX);
        SIMULATION_DATA.SNSTA_1(simulation_step,:) = [simulation_time; x_asnsta_1; control_asnsta_1; disturbance; epsilon_asnsta_1; d_est_x_asnsta_1];
        dot_x_asnsta_1 = control_asnsta_1 + disturbance;
        x_asnsta_1 = x_asnsta_1 + dot_x_asnsta_1*PARAMETERS.SIMULATION.SAMPLING_TIME;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------------------------------------------------------------------------------------------------% 
        %----------------------------------------------------- SNSTA_2 --------------------------------------------------------% 
        %---------------------------------------------------------------------------------------------------------------------% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PARAMETERS.SIMULATION.NOISE_MODULE_DB > 0
            noise_x_asnsta_2 = awgn(x_asnsta_2 , PARAMETERS.SIMULATION.NOISE_MODULE_DB, "measured");   
        else
            noise_x_asnsta_2 = x_asnsta_2;
        end
        z_asnsta_2 = noise_x_asnsta_2 - est_x_asnsta_2;
        dot_est_x_asnsta_2 = lambda*z_asnsta_2 + beta*tanh(gamma*z_asnsta_2) + epsilon_asnsta_2 + control_asnsta_2;
        dot_epsilon_asnsta_2 = ((lambda^2)/4.0)*z_asnsta_2 + kappa*tanh(gamma*z_asnsta_2);
        epsilon_asnsta_2 = epsilon_asnsta_2 + dot_epsilon_asnsta_2*PARAMETERS.SIMULATION.SAMPLING_TIME;
        est_x_asnsta_2 = est_x_asnsta_2 + dot_est_x_asnsta_2*PARAMETERS.SIMULATION.SAMPLING_TIME;
        d_est_x_asnsta_2 = lambda*z_asnsta_2 + beta*tanh(gamma*z_asnsta_2) + epsilon_asnsta_2;
        U = PARAMETERS.CONTROL.MAX - abs(d_est_x_asnsta_2);
        control_asnsta_2 = -U*tanh(PARAMETERS.CONTROL.CHI_2*noise_x_asnsta_2) - d_est_x_asnsta_2;
        control_asnsta_2 = function_sat(control_asnsta_2, PARAMETERS.CONTROL.MAX);
        SIMULATION_DATA.SNSTA_2(simulation_step,:) = [simulation_time; x_asnsta_2; control_asnsta_2; disturbance; epsilon_asnsta_2; d_est_x_asnsta_2];
        dot_x_asnsta_2 = control_asnsta_2 + disturbance;
        x_asnsta_2 = x_asnsta_2 + dot_x_asnsta_2*PARAMETERS.SIMULATION.SAMPLING_TIME;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update time
        simulation_time = simulation_time + PARAMETERS.SIMULATION.SAMPLING_TIME;
    end
end