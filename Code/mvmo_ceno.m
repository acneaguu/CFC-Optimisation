%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%By: Dr. Jose L. Rueda
%%Date: 02.03.2012
%%Edited: Alex Neagu & Jinhan Bai
%%Date: 11/05/2020

%%This is the source code for the MVMO-SHM algorithm. It will return the
%%solution resulting in the best fitness of the function fhd. The
%%boundaries of the optimisation variables are given by lb and ub. Before
%%running this function, initialise using initialise_mvmoshm function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Fout,Xout] = mvmo_ceno(fhd,lb,ub)
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	global proc
    global ps 
    global parameter 
    global bests
    global fitness
    global changed_best
    global x_normalized_save
    global merke 
    global shape
    global meann
    global DD
    global ipp
    global updated
    global nnnnnn
    
    % put the bounds into the corresponding internal variables
    ps.x_min = lb;
    ps.x_max = ub;
    
    % Reinitialization of local
    % function evaluation counter.
    proc.i_eval=0;

    % Function evaluation at which
    % the last update of the global
    % best solution occurred.
    % Refers to the internal evaluation
    % using static penalty constraint 
    % handling method.
	proc.last_improvement=1;
    % Signalizing your 
    % to stop running.
    proc.finish=0;   
    % Dimensionality of test case.
    D=ps.D;
    DD=1:D;
    
    
    n_eval =proc.n_eval;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    min_eval_LS =round(parameter.ratio_local*proc.n_eval);           %  minimum number of runs without local search
    max_eval_LS =round(1.00*proc.n_eval+10);
    derivative='forward';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    alpha=1d0    ;                        
    n_steps=1   ;
    bestp=1;
    local_i=0;
    parameter.scaling=(ps.x_max-ps.x_min);
    %%----------------- Create initial random population ------------------   
    merke=zeros(parameter.n_par,ps.D);
    xx=zeros(parameter.n_par,ps.D);
    x_norm=xx;
    for iijj=1:parameter.n_par %Start initialization: x_normalized
          for jjkk=DD
              xx(iijj,jjkk)=ps.x_min (jjkk) + rand(1)*(ps.x_max(jjkk)-ps.x_min(jjkk));
          end
        x_norm(iijj,:)=(xx(iijj,:)-ps.x_min)./parameter.scaling;
    end % End initialization
    x_normalized=x_norm;
    x_normalized_save=x_norm;
    
    %% ------------------ Initialize control parameters --------------------
    n_par = parameter.n_par; 
    n_to_save = parameter.n_tosave; 
    fs_factor_start = parameter.fs_factor_start;
    fs_factor_end = parameter.fs_factor_end;
    IX=ones(n_par); 
    shape = zeros(n_par,D);
    shape(:,:) =1.d0;
    for i=1:n_par
     IX(i)=i;
    end
%     bestp=IX(1);
    ratio_gute_min=parameter.ratio_gute_min;
    ratio_gute_max=parameter.ratio_gute_max;
    n_randomly_ini = parameter.n_random_ini;
    n_randomly_last = parameter.n_random_last;
    l_vari=n_to_save  ; % can be changed
      if  l_vari<2 || l_vari > n_to_save 
        l_vari=n_to_save ;
      end  
    
    %% --------------------- Data structure for the table ---------------------
    bests = NaN*zeros(parameter.n_tosave,ps.D,n_par,'double');
    fitness = Inf*ones(parameter.n_tosave,n_par,'double');
    
    %% ----------------------------- Mapping ----------------------------------
    updated=zeros(n_par,D);
    nnnnnn=zeros(n_par,1);
    values_noneq=zeros(parameter.n_tosave,1);
    meann = x_normalized;
    meann(:,:)= 0.5d0  ;
    local_search=zeros(n_par,1);
    considered = zeros(1,ps.D);
    
   %% ---------------------- Checking correctness ------------------------   
   
    if (n_randomly_last<1)
        n_randomly_last=1;
    end
    
    if (n_randomly_last>ps.D)
        n_randomly_last=ps.D;
    end
    
    if (n_randomly_ini>ps.D)
        n_randomly_ini=ps.D;
    end  
    if (n_eval<=0)
        n_eval=100000d0;
    end

    if (fs_factor_start<=0.0)
        fs_factor_start=1d0;
    end
    
    if (fs_factor_end<=0)
        fs_factor_end=1d0;
    end

    fs_factor0=fs_factor_start;
    
    yes_n_randomly=true;
    if (n_randomly_ini==n_randomly_last)
        n_randomly_n=n_randomly_last;
        yes_n_randomly=false;
    end
    
    if (n_to_save<=1)
        n_to_save=2d0;
    end
    
    yes_fs_factor=true;
    if (fs_factor_start==fs_factor_end)
        yes_fs_factor=false;
    end

%% ----------------------------- Counters ---------------------------------

    no_in = zeros(1,n_par);
    no_inin = zeros(1,n_par);
    goodbad=ones(n_par);
    firsttime=true;
    f_LocAlgo=1;
     min_eval_LS =round(min_eval_LS/f_LocAlgo);

    delta_nrandomly=n_randomly_ini-n_randomly_last ; 
  while 1       
        %Evaluating the particles.....
            ff=proc.i_eval/n_eval   ;
            ff2=ff*ff;
            ff4=ff2*ff2;
            eins_minus_ff2=1.d0-ff2;
            eins_minus_ff2_q=eins_minus_ff2*eins_minus_ff2;
            eins_minus_ff2_q_p=eins_minus_ff2_q+1.d-4;
            vvqq=10.d0^(-(3.5d0+ff*5.d0));
            %Determining the relative number of particles belonging to the group of
            %good particles
            border_gute0=ratio_gute_max-ff*(ratio_gute_max-ratio_gute_min);
            border_gute=round((n_par*border_gute0)) ;
             if border_gute < 3 && n_par > 3
                 border_gute=3;
             end
            if border_gute > n_par
                border_gute=n_par ;
             end                     
            %Selecting the subset of variables to be mutated
            if yes_n_randomly
                n_randomly_X=round(n_randomly_ini-ff4*delta_nrandomly) ;
                n_randomly_n=round(n_randomly_last+rand(1)*(n_randomly_X-n_randomly_last)) ;
                if rand(1)>0.75
                    n_randomly_n=round(n_randomly_n/2);
                end
            end
            %Quadratic variation of fs_factor0
            if yes_fs_factor
                fs_factor0=fs_factor_start+ff*(fs_factor_end-fs_factor_start)  ;
            end
            
      
     ipx=0 ;
     while ipx  < n_par 
           ipx=ipx+1 ;
           ipp=IX(ipx);
% -------------MAPPING----------------
%   only for not local search runs
       considered(:)=0;
       if nnnnnn(ipp) >=n_to_save    
          if n_randomly_n>0 && n_randomly_n<=D 
                D_imc=ps.D;%-i_not;
               getestet=0;
               in_randomly=0;
               while in_randomly<n_randomly_n &&  getestet<D_imc
                izmm=randi(ps.D); %round(rand(1)*(n_var-1))+1;  sort
                izzz=izmm;
                isrepeat = true;
                in_randomlyx=0;
                while   isrepeat &&  getestet<D_imc  
                    if (izzz< 1) 
                         izzz=ps.D;
                    end
                    if  considered(izzz)==0  
                    in_randomlyx=in_randomlyx+1;
                    in_randomly=in_randomly+1;
                    getestet=getestet+1;
                    considered(izzz)=1;
                    izzz=izzz-1;
                    if in_randomlyx < n_steps  && in_randomly < n_randomly_n  
                        isrepeat = true;
                    else
                        isrepeat = false;
                    end
                    else
                     izzz=izzz-1;   
                   end
                end
                end
          end
%
            bestp=IX(1);
            if ipx==1
            onep1=IX(1);
            worstp=IX(1);
            elseif ipx==2 
            onep1=IX(2);
            worstp=IX(2);
            elseif ipx==3 
            onep1=IX(2);
            worstp=IX(3);
            else
            onep1=IX(randi(border_gute-2)+1);
            worstp=IX(border_gute);
            end
            bbb=1.1d0+(rand(1)-0.5d0)*2.0d0  ; 
            beta1 =alpha*3.d0*bbb*((1.d0+0.0*ff2)*rand(1) - (1.d0-ff2)*0.9d0)*rand(1)  ;  % the best one
            gut=false;
            if goodbad(ipp)  %&& rand(1)   
                gut=true;
            end
         x_normalized(ipp,1:D)= bests(1,1:D,bestp) ; 
         for ivar=DD
             if considered(ivar) == -1  
                 continue
             end
               if considered(ivar) == 1   
                     if rand(1) > ff    
                      xxxx = rand(1)  ;    
                    else
                      xxxx = -1.  ; 
                       while (xxxx > 1.d0 || xxxx < 0.d0) 
%                             xxxx=normalfunc(0.5d0 , eins_minus_ff2_q_p )  ;     
                             xxxx =0.5d0+eins_minus_ff2_q_p*(sum(rand(1,12))-6.d0)/6.d0;
                       end  
                    end  
                       fs_factor1=abs(fs_factor0*(4.d0 + 1.65d0*(rand(1)-0.15))) ; 
                       if fs_factor1 < ff
                         fs_factor1=ff;
                       end
% Update the required mean   'meann(izuff,ivar' and 'shape(izuff,ivar)'
     izuff=1 ;
     if border_gute>1 izuff=IX(randi(border_gute)); end;
     for i_update=[izuff,ipp]
       if updated(i_update,ivar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            iz=1;
            wert=bests(1,ivar,i_update);
            values_noneq(iz )=wert;
%             vmean=wert;
            for ii_jj=2:nnnnnn(i_update)
                izz = iz;
                thesame = false;
                for kk_ii=1:izz
                    wert=bests(ii_jj,ivar,i_update);
                    if  abs(values_noneq(kk_ii ) - wert) < vvqq  
                        thesame=true;
%                         values_noneq(kk_ii )=wert;
                        break;
                    end
                end
                if ~thesame 
                    iz = iz+1;
                    values_noneq(iz )=wert;
                end
            end
                    vmean=sum(values_noneq(1:iz))  ; 
            if (iz>1)
                   vmean=vmean/iz;
                   meann(i_update,ivar)= vmean;
                   values_noneq(1:iz)=values_noneq(1:iz )-vmean; 
                   vv =sum(values_noneq(1:iz ).*values_noneq(1:iz ))/iz;
                   if vv > 1.d-50 
                      shape(i_update,ivar)=-reallog(vv);
                   end
                   updated(i_update,ivar)=false;
            end
       end
%        end
     end
%  H-function calculation included in the main code
       if  xxxx <= 0.5d0 
       %  %Forward
       sss1 =  shape(ipp,ivar)*fs_factor1;
       Mittelwert=meann(izuff,ivar);
       if n_par==1   &&  rand(1) > 0.5d0
           Mittelwert=bests(randi(nnnnnn(izuff)),ivar,izuff);%Mittelwert+(0.5d0-rand(1))*2.55d0;
           if Mittelwert > 1.0d0
               Mittelwert=0.5d0;
           elseif Mittelwert < 0.0d0
               Mittelwert=0.5d0;
           end
       end
%             
       zw1=1.d0-Mittelwert; 
          if  zw1 > 1.d-5 
             s11=sss1/zw1 ;
          else
              Mittelwert=0.99999d0; 
              zw1=1.d0-Mittelwert; 
              s11=sss1/zw1 ;
          end
         hmitte_f=-Mittelwert/(0.5d0*s11+1.d0)+Mittelwert ;
         hforw=-Mittelwert/(xxxx*s11+1.d0)+Mittelwert;
         hcorr=(Mittelwert-hmitte_f)*(xxxx+xxxx);
         xxxx=hforw+hcorr ;
       else 
     %    %Backward
         sss1 =  shape(izuff,ivar)*fs_factor1;
         Mittelwert=meann(izuff,ivar);
       if n_par==1  &&  rand(1) > 0.5d0
%            Mittelwert=bests(randi(nnnnnn(izuff)),ivar,izuff);%Mittelwert+(0.5d0-rand(1))*2.55d0;
           Mittelwert=bests(randi(nnnnnn(izuff)),ivar,izuff);%Mittelwert+(0.5d0-rand(1))*2.55d0;
           if Mittelwert > 1.0d0
               Mittelwert=0.5d0;
           elseif Mittelwert < 0.0d0
               Mittelwert=0.5d0;
           end
       end
         if  abs(Mittelwert) > 1.d-5   
           zw1=1.d0-Mittelwert; 
           s11=sss1/Mittelwert ;
         else
           Mittelwert=0.00001d0 ;
           zw1=1.d0-Mittelwert; 
           s11=sss1/Mittelwert ;
         end
         zw2=1.d0-xxxx;
         hmitte_b=zw1/(0.5d0*s11+1.d0);  
         hback=zw1/(zw2*s11+1.d0) + Mittelwert;
         hcorr=hmitte_b*(zw2+zw2);
         xxxx=hback-hcorr; 
       end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
               else
                        if gut
                             if n_par > 1
                              xxxx =0.99999*bests(1,ivar,ipp)+ rand(1)*0.00001*bests(1,ivar,bestp);  
                             else
                              xxxx =bests(1,ivar,ipp) ;  
                             end
                             
                        else
                            if rand(1)>ff
                              xxxx =bests(1,ivar,onep1)+ (beta1)*(bests(1,ivar,bestp)-bests(1,ivar,worstp));   
                            else
                              if bestp ==  onep1
                                aaaa=1.d0-ff2;
                                bbbb=(0.5d0-rand(1))*aaaa;
                                xxxx =bests(1,ivar,onep1)+ bbbb*bests(1,ivar,onep1);   
                              else
                               aaaa=rand(1)*0.01d0;
                               bbbb=1.d0-aaaa;
                               xxxx =bbbb*bests(1,ivar,onep1)+ aaaa*bests(1,ivar,bestp);   
                              end
                            end
                        end
               end
                       if xxxx > 1.0d0    
                          xxxx=1.0d0;
                       elseif xxxx < 0.0d0  
                          xxxx=0.0d0 ;
                       end
                       x_normalized(ipp,ivar)=xxxx;            
         end
       else
          x_normalized(ipp,:)=rand(ps.D,1)  ;%table.bests(1,DD,ipp); %Local best-based parent assignment during independent evaluations
       end
% -------------END of MAPPING----------------
%Denormalize to the original [min, max] range 
             if rand < parameter.local_prob &&  proc.i_eval > min_eval_LS &&  proc.i_eval < max_eval_LS && ipp==bestp 
                      if local_i==0
                          x_normalized_save(ipp,:)=bests(1,:,bestp) ;
                          ffbest=fitness(1:1,bestp);
                          if ffbest < 1.d9  
                              optimmethod='active-set';
                          else 
                              optimmethod='interior-point';
                          end   
                              
                      else
                         x_normalized_save(ipp,:)= x_normalized(ipp,:);
                      end
                     x_normalized_save(ipp,:) = ps.x_min+parameter.scaling.* x_normalized_save(ipp,:);
                 [ffx,~,~,x_normalized(ipp,:),FEVALS] = LocalSearchMVMOSH(fhd,x_normalized_save(ipp,:),proc.n_eval -proc.i_eval); %Local search
                 local_search(ipp)=local_search(ipp)+1; %COUNTING HOW MANY TIMES LOCAL SEARCH WAS CALLED
                 x_normalized(ipp,:) = (x_normalized(ipp,:)-ps.x_min)./parameter.scaling;
                  local_i=local_i+1  ;
                  n_randomly_ini=n_randomly_last;

             else
                 x_normalized(ipp,:) = ps.x_min+parameter.scaling.* x_normalized(ipp,:);
                 [ffx,~,~,~]=feval(fhd,x_normalized (ipp,:)); %Problem evaluation
                 x_normalized(ipp,:) = (x_normalized(ipp,:)-ps.x_min)./parameter.scaling;
            end
                 fitness_new=ffx;%oox+1.d-12*(ffx-oox); % Constraint handling outsourced so far   
            
            if proc.finish
                Fout = min(fitness(1,:));
                bp_index = fitness(1,:)==Fout;
                Xout = bests(1,:,bp_index).*parameter.scaling +ps.x_min;
                if size(Xout,3) > 1
                   Xout = Xout(:,:,1); 
                end
                return;
            end
%
            Fill_solution_archive(fitness_new) ;
           
            %Determining the proportion of good particles
                   if n_par > 1
                    if  firsttime
%                         A(1:n_par)=fitness(1,1:n_par);
                      firsttime=false  ;
                      [~,IX] = sort(fitness(1,1:n_par));
                      goodbad(1:n_par)=0;
                      goodbad(IX(1:border_gute))=1;
                    
                    elseif changed_best
                        [~,IX] = sort(fitness(1,1:n_par));
                      goodbad(1:n_par)=0;
                      goodbad(IX(1:border_gute))=1;
                    end
%                     end
                   end
                     
                   ipp=IX(ipx); 
                   if n_par < 4 
                    goodbad(ipp)=1   ;
                   end
%                 end
     end %End n_par loop
end %End while loop        

     %% ----------------------- Complementary functions ------------------------
     
%         function [lnnnnnn,lupdated]=Fill_solution_archive (lnnnnnn,lupdated)
        function Fill_solution_archive(fitness_new)
%
        no_in(ipp) = no_in(ipp)+1;
        changed = false;
        changed_best=false;
        if no_in(ipp) ==1 % the solution coming to the table for the first time large
            fitness(1:n_to_save,ipp) = 1.e200;
            bests(1,:,ipp)   = x_normalized(ipp,:) ;
            fitness(1,ipp) = fitness_new ;
            no_inin(ipp)=no_inin(ipp)+1;
            changed_best=true;
            nnnnnn(ipp) = 1;
%             if no_inin(ipp) >= n_to_save
%                 updated(ipp,DD)=true;
%             end
        else % not for the first time and check for the update
           i_position =0;
               for ij=1:n_to_save 
                   if fitness_new < fitness(ij,ipp) %&& xt.feasibility == table.feasibility(ij,:,ipp)) || (table.feasibility(ij,:,ipp) <  xt.feasibility)                                 
                       i_position = ij;
                       changed =true;
                         if (ij<n_to_save)
                           no_inin(ipp) = no_inin(ipp)+1;  % how many times good solutions were found   
                           if no_inin(ipp) >= n_to_save
                             updated(ipp,:)=true;
                           end
%                          else
%                              dsdsd=0;
                         end
                       break;
                   end
               end
        end

        if changed   % if the new individual is better than any archived individual.
                     % Move the individuals and corresponding fitness values
                     % downward so that the individuals are sorted based on the
                     % fitness value in a descending order             
            if i_position==1
              changed_best=true;
            end     
              i_end=nnnnnn(ipp)+1; %i_position+1;
              if i_end > n_to_save
                  i_end =n_to_save;
              end
              if   i_end~=i_position
              for i_lauf=i_end:-1:i_position+1
                bests(i_lauf,:,ipp) = bests(i_lauf-1,:,ipp);
                fitness(i_lauf,ipp)= fitness(i_lauf-1,ipp);
              end
            % save the new best
              bests(i_position,:,ipp) = x_normalized(ipp,:);
              fitness(i_position,ipp) = fitness_new;
              
              else
        
                bests(i_position,:,ipp) = x_normalized(ipp,:);
                fitness(i_position,ipp) = fitness_new;
              end
              if nnnnnn(ipp) < n_to_save 
                  if fitness(nnnnnn(ipp)+1,ipp) < 1.e200
                    nnnnnn(ipp)=nnnnnn(ipp)+1;
                  end
              end
        end
        end
    

         %% ----------------------- Complementary functions ------------------------
     function [ffx,oox,ggx,xn_out,FEVALS] = LocalSearchMVMOSH(fhd,xx_yy,FEsAllowed)
        global PPL GGL
        if FEsAllowed <= 0, return, end
        options = optimoptions(@fmincon,'Display','none','algorithm',optimmethod,'UseParallel','never','MaxFunEvals',FEsAllowed,'FinDiffType',...
                          derivative,'TolFun',1.d-3,'TolX',1d-10,'HessianApproximation','lbfgs','DiffMinChange',1.d-10);
        [Xsqp, FUN , ~ , output]=...
            fmincon(@(xx_yy)LSearch(xx_yy,fhd),xx_yy,[],[],[],[],ps.x_min,ps.x_max,[],options);

        
        FEVALS=output.funcCount  ;
        for nvar=1:size(xx_yy,2)
            if isnan(Xsqp(1,nvar))
                Xsqp=xx_yy;
                break;
            end
        end
        
        xn_out = Xsqp;
        ffx=FUN;
        oox=PPL; 
        ggx=GGL;
    end

    function J=LSearch(xx_yy2,fhd)
        global PPL GGL 
        [J,PPL,GGL,~] = feval(fhd,xx_yy2);
        
    end
        

    
    
    
end