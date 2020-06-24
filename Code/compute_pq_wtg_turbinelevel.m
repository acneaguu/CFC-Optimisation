function [P,Q] = compute_pq_wtg_turbinelevel(windspeed)
    %%wtg data
    v_c_in = 3 ; %general cut in speed in m/s
    v_r = 14 ; %general rated speed in m/s, fixed for step 2a
    v_c_off = 25 ; %general cut off speed in m/s
    P_wt_max = [4.2 4 4.2 2.35]; %Rated power in MW for each string
    Q_wt_max = [2.65 2.8 2.8 0.85] ; %Max positive/negative reactive power in MVAr generated for eacht string 1-13
    
    %% Compute P
    %%init
    nsamples = length(windspeed);
    P = zeros(nsamples, 4); %Make matrix for generated nominal wind speed

    %%convert windspeed to Power using cubic wind power model
    for i=1:4 %makes a 1x13
        for j=1: nsamples %makes a 12x1
        %%apply boundary conditions to determine which equation is valid
            if(windspeed(j) <= v_c_in) 
                P(j,i) =0; 
            elseif ( (windspeed(j) > v_c_in) &&  (windspeed(j) <= v_r) )    
                P(j,i)=P_wt_max(i)*(windspeed(j)^3-v_c_in^3)/(v_r^3-v_c_in^3);    
            elseif ((windspeed(j) > v_r) && (windspeed(j) <= v_c_off))
                P(j,i) = P_wt_max(i); %assigns every power to every dispatch     
            elseif (windspeed(j) > v_c_off)
                P(j,i) = 0;
            end 
        end
    end
    
    
    
    %% Compute Q
    
%     rc_string_in = [6.625  6.625   11.129  13.625  6.22    11.129  6.22    6.22    11.2    7.36    6.22    12.25 10.64]; %specifies slope MVAr/MW at beginning for each WTG string
%     P_reg_in = [0.1 0.1 0.0544  0.047619    0.1 0.0544  0.1 0.1 0.0606  0.061   0.1 0.0556  0.0556]; %Percentage of total power to reach Q max in capability curve (per turbine)
%     rc_string_end = -[1.5    1.5     0.507   3.143   1.48   0.507   1.8519  1.48    0.4408  0.4047   1.48    0.531   0.531]; %specifies slope MVAr/MW at end for each WTG string
%     %new_Q_max = [0.72   0.72    0.5117   0.512   0.65    0.5117   0.59    0.643   0.547   0.302   0.643 0.4898   0.4898]' * Q_wt_max; %Calculates new available MVAr at P_wt_max
%     P_reg_end = [0.88   0.88    0.372   0.917   0.83125   0.372   0.83125   0.83125   0.3023 0.225  0.83125   0.346   0.346]; %Percentage of total power to reach final Q in capability curve

rc_string_in = [6.625 6.2222 7 4.25];
rc_string_end = -[2.095 1.7037 1.357 4.3];
P_reg_in = [(0.4/4.2) (0.4/4) (0.4/4.2) (0.2/2.35)];
P_reg_end = [(3.675/4.2) (3.325/4) (3.85/4.2) (2.25/2.35)];

    Q = zeros (nsamples, 4);

    for i=1:4
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