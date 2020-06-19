%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Original authors: Farley Rimon & Marouane Mastouri
%%This is a modified version.

%%README:
%%This function is used to compute the active and reactive power outputs of
%%the different strings as function of the wind speed. The output are
%%two vectors containing the P and Q of each string. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,Q] = compute_pq_wtg(windspeed)
    %%WTG data
    %General cut in speed in m/s
    v_c_in = 3;    
    
    %General rated speed in m/s
    v_r = 14; 
    
    %General cut off speed in m/s
    v_c_off = 25;  
    
    %Rated power in MW for each string
    P_wt_max = [33.6 29.4 29.4 16.8 32 29.4 16 28 33 29.15 28 28.8 28.8]; 
    
    %Max positive/negative reactive power in MVAr generated for each string
    Q_wt_max = [21.2 18.55 19.15 10.9 22.4 19.15 12.2 19.6 22.4 13.248 19.6 19.6 19.6] ; 
    
    %% Compute P
    %%Initialise P
    nsamples = length(windspeed);
    P = zeros(nsamples,13);

    %%Convert windspeed to Power using cubic wind power model
    for i=1:13 
        for j=1: nsamples
            %%Apply boundary conditions to determine which equation is valid
            if(windspeed(j) <= v_c_in) 
                P(j,i) =0; 
            elseif ( (windspeed(j) > v_c_in) &&  (windspeed(j) <= v_r) )    
                P(j,i)=P_wt_max(i)*(windspeed(j)^3-v_c_in^3)/(v_r^3-v_c_in^3);    
            elseif ((windspeed(j) > v_r) && (windspeed(j) <= v_c_off))
                P(j,i) = P_wt_max(i);    
            elseif (windspeed(j) > v_c_off)
                P(j,i) = 0;
            end 
        end
    end
    %% Compute Q
    %%Initialise Q
    Q = zeros (nsamples, 13);
    
    %%Specifies slope MVAr/MW at beginning for each WTG string
    rc_string_in = [6.625  6.625   11.129  13.625  6.22    11.129  6.22    6.22    11.2    7.36    6.22    12.25 10.64]; 
    
    %%Percentage of total power to reach Q max in capability curve (per turbine)
    P_reg_in = [0.1 0.1 0.0544  0.047619    0.1 0.0544  0.1 0.1 0.0606  0.061   0.1 0.0556  0.0556]; 
    
    %%Specifies slope MVAr/MW at end for each WTG string
    rc_string_end = -[1.5    1.5     0.507   3.143   1.48   0.507   1.8519  1.48    0.4408  0.4047   1.48    0.531   0.531]; 
    
    %%Percentage of total power to reach final Q in capability curve
    P_reg_end = [0.88   0.88    0.372   0.917   0.83125   0.372   0.83125   0.83125   0.3023 0.225  0.83125   0.346   0.346]; 

    %%Compute Q depending on the region in the  capability curve
    for i=1:13
        for j=1:nsamples  
            if (P(j,i) < (P_reg_in(i)*P_wt_max(i)))
                Q(j,i) = rc_string_in(i)*P(j,i);
                if (Q(j,i) >= Q_wt_max(i))
                    Q(j,i) = Q_wt_max(i);
                end 
            elseif ((P(j,i) >= (P_reg_in(i)*P_wt_max(i))) && (P(j,i) < (P_reg_end(i)*P_wt_max(i))))
                Q(j,i) = Q_wt_max(i);
            elseif ((P(j,i) >= (P_reg_end(i)*P_wt_max(i))))
                Q(j,i) = Q_wt_max(i) + rc_string_end(i).*(P(j,i)-P_reg_end(i)*P_wt_max(i));
                if (Q(j,i) < 0)
                    Q(j,i) = 0;
                end 
            end 
        end
    end
end