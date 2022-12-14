Module mcmc_uncorr

contains

SUBROUTINE hbtdmcmc(beta,k, x, h_agg,y,n0, plon, plat, regions_v, &
pdf_param1, pdf_param2, lons,lats, h_v, sigma_model, sigma_measure, &
R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, &
sigma_clon, sigma_clat, rjmcmc, para_temp, & 
lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in, &
pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, pdf_param2_pdf, &
stepsize, stepsize_pdf_p1,stepsize_pdf_p2, nIt, nsub, nit_sub, nIC, &
nbeta, kmax, kICmax, nmeasure, Ngrid, ydim1, ydim2, nIC1, &
k_out, x_out, y_out, regions_out, plon_out, plat_out, sigma_model_out, sigma_y_out, &
n0T_out, pdf_param1_out, pdf_param2_out, accept, reject, &
accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, &
accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, &
tot_acc_x, tot_acc_p1, tot_acc_p2, tot_acc_sigma_y, &
accept_all, reject_all, accept_birth_all, reject_birth_all, &
accept_death_all, reject_death_all, accept_move_all, reject_move_all, &
accept_sigma_y_all, reject_sigma_y_all)


IMPLICIT NONE

! THIS SUBROUTINE IS TO PERFORM A TRANSDIMENSIONAL INVERSION USING REVERSIBLE JUMP MCMC
! THIS CODE RELIES ON MEASUREMENTS BEING UNCORELLATED
! THIS IS A DANGEROUS ASSUMPTION IF APPLIED POORLY BUT ALLOWS MUCH FASTER CALCULATION
! THAT MEANS THIS SCRIPT IS GOOD FOR TESTING STUFF 

! REQUIRES:
! Inputs to be genrated in Python using script tdmcmc_template.py

! OUTPUTS:
! Will be passed back to tdmcmc_template.py

! PARALLEL TEMPERING:
! The inversion can be performed with parallel tempering or without
! The default is without, performed on a single processor

! The code is set up such that Voronoi nuclei can only lie on the underlying grid
! i.e. Only ever on the centre points of each grid cell

!! INPUTS !!!!!!!!!!!!!!!!!!!
! Dimensions
INTEGER nbeta
INTEGER kmax
INTEGER kICmax
INTEGER nmeasure
INTEGER Ngrid
INTEGER nit_sub
INTEGER nIt
INTEGER burn_in
INTEGER nsub
INTEGER kmin
INTEGER nIC
INTEGER ydim1
INTEGER ydim2
INTEGER nIC1
! Single Variables
REAL lonmin
REAL lonmax
REAL latmin
REAL latmax
REAL sigma_bd
REAL sigma_clon
REAL sigma_clat
INTEGER sigma_model_pdf
INTEGER pdf_param1_pdf
INTEGER pdf_param2_pdf
INTEGER rjmcmc
INTEGER para_temp
! Input arrays
REAL stepsize(nIC1)
REAL stepsize_pdf_p1(nIC1)
REAL stepsize_pdf_p2(nIC1)
INTEGER x_pdf_all(nIC1)
REAL beta(nbeta)
INTEGER k(nbeta)
REAL x(kICmax, nbeta)               
REAL pdf_param1(kICmax,nbeta)            
REAL pdf_param2(kICmax,nbeta)                
REAL pdf_p1_hparam1(nIC1)
REAL pdf_p1_hparam2(nIC1)
REAL pdf_p2_hparam1(nIC1)
REAL pdf_p2_hparam2(nIC1)
REAL h_agg(nmeasure,kICmax,nbeta)
REAL y(nmeasure) 
REAL n0(nmeasure, nbeta) 
REAL n0T(nbeta)
REAL sigma_y(nmeasure,nbeta)
REAL sigma_model(ydim2, nbeta)
REAL stepsize_sigma_y(ydim2)
INTEGER R_indices(ydim1,ydim2)
REAL sigma_measure(nmeasure)
REAl sigma_model_hparam1(ydim2)
REAl sigma_model_hparam2(ydim2)
REAL plon(kmax,nbeta)
REAL plat(kmax,nbeta)
INTEGER regions_v(Ngrid,nbeta)
REAL lons(Ngrid)
REAL lats(Ngrid)
REAL h_v(nmeasure, Ngrid)
! Outputs
INTEGER k_out(nit_sub)
REAL x_out(kICmax,nit_sub)  
REAL y_out(nmeasure,nit_sub)  
REAL pdf_param1_out(kICmax,nit_sub)             
REAL pdf_param2_out(kICmax,nit_sub)   
REAL plon_out(kmax,nit_sub)
REAL plat_out(kmax,nit_sub)
INTEGER regions_out(Ngrid,nit_sub)
REAL sigma_y_out(nmeasure, nit_sub)
REAL sigma_model_out(ydim2,nit_sub)
REAL n0T_out(nit_sub)
REAL tot_acc_x(nIC1)
REAL tot_acc_p1(nIC1)
REAL tot_acc_p2(nIC1)
REAL tot_acc_sigma_y(ydim2)
INTEGER accept(nIC1)
INTEGER reject(nIC1)
INTEGER accept_birth, reject_birth
INTEGER accept_death, reject_death, accept_move, reject_move
INTEGER accept_swap, reject_swap
INTEGER accept_sigma_y, reject_sigma_y
INTEGER accept_all(nIC1,nbeta)
INTEGER reject_all(nIC1,nbeta)
INTEGER accept_birth_all(nbeta), reject_birth_all(nbeta)
INTEGER accept_death_all(nbeta), reject_death_all(nbeta)
INTEGER accept_move_all(nbeta), reject_move_all(nbeta)
INTEGER accept_sigma_y_all(nbeta), reject_sigma_y_all(nbeta)
! INTERMEDIATE VARIABLES
INTEGER it, ibeta, remain_it, pair1,pair2, ib, it_sub, remain, kIC       !remain_dim
INTEGER remain_swap,jj
REAL u, u1,u2, randomu,pT_chain, beta1,beta2
INTEGER k_it(nit_sub)
REAL x_it(kICmax,nit_sub)    
REAL y_it(nmeasure,nit_sub)                           
REAL pdf_param1_it(kICmax,nit_sub)   
REAL pdf_param2_it(kICmax,nit_sub)   
REAL plon_it(kmax,nit_sub)               
REAL plat_it(kmax,nit_sub)              
REAL sigma_model_ap(ydim2)           
INTEGER regions_it(Ngrid,nit_sub)
REAL sigma_y_it(nmeasure,nit_sub), sigma_model_it(ydim2,nit_sub)     
REAL n0T_it(nit_sub)
REAL sigma_y_temp(nmeasure), y_error_temp(nmeasure)
INTEGER acc_h_batch(nIC1)
INTEGER rej_h_batch(nIC1)
INTEGER accept_batch(nIC1)
INTEGER reject_batch(nIC1)
INTEGER acc_y_batch(ydim2)
INTEGER rej_y_batch(ydim2)
! SUBROUTINE INPUTS
REAL betaib, n0Tib
REAL xib(kICmax), plonib(kmax), platib(kmax), n0ib(nmeasure)           
REAL pdf_param1ib(kICmax), pdf_param2ib(kICmax)
REAL h_aggib(nmeasure,kICmax)
REAL sigma_yib(nmeasure), sigma_modelib(ydim2)
INTEGER kib
INTEGER regions_vib(Ngrid)
! SUBROUTINE OUTPUTS
REAL n0Tib1                     
REAL xib1(kICmax), plonib1(kmax), platib1(kmax), n0ib1(nmeasure)      
REAL pdf_param1ib1(kICmax), pdf_param2ib1(kICmax)
REAL h_aggib1(nmeasure,kICmax)
REAL sigma_yib1(nmeasure), sigma_modelib1(ydim2)
INTEGER kib1, rejectib1, acceptib1, reject_yib1, accept_yib1
INTEGER acceptxib1(nIC1), rejectxib1(nIC1), acc_bxib1(nIC1), rej_bxib1(nIC1)
INTEGER acc_byib1(ydim2), rej_byib1(ydim2)
INTEGER regions_vib1(Ngrid)
REAL detval(nbeta)
REAL detvalib, detvalib1
REAL stepsize_sig_ib1(ydim2)
REAL stepsize_ib1(nIC1)     
REAL stepsize_p1_ib1(nIC1)
INTEGER acc_prob_p1_ib1(nIC1)
REAL stepsize_p2_ib1(nIC1)
INTEGER rej_prob_p1_ib1(nIC1)

!! F2PY IN/OUT COMMANDS !!!!!!!
!f2py intent(in) beta,k, x, h_agg,y,n0, plon, plat, regions_v
!f2py intent(in) pdf_param1, pdf_param2, lons,lats, h_v, sigma_clon, sigma_clat  
!f2py intent(in) sigma_model, sigma_measure, rjmcmc, para_temp
!f2py intent(in) R_indices, sigma_model_hparam1, sigma_model_hparam2
!f2py intent(in) stepsize_sigma_y, sigma_model_pdf
!f2py intent(in) lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in
!f2py intent(in) pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf
!f2py intent(in) pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf
!f2py intent(in) stepsize, nIt,nsub,nit_sub, nIC1
!f2py intent(in) nIC, nbeta, kmax, kICmax, nmeasure, Ngrid, ydim1, ydim2
!f2py intent(out) k_out, x_out, y_out, regions_out, plon_out, plat_out
!f2py intent(out) sigma_model_out, sigma_y_out
!f2py intent(out) n0T_out, pdf_param2_out, pdf_param1_out, reject_swap
!f2py intent(out) accept, reject, accept_swap, accept_birth, accept_death, accept_move
!f2py intent(out) reject_birth, reject_death, reject_move, accept_sigma_y, reject_sigma_y
!f2py intent(out) tot_acc_sigma_y, tot_acc_x, tot_acc_p1, tot_acc_p2
!f2py intent(out) accept_all, reject_all, accept_birth_all, reject_birth_all
!f2py intent(out) accept_death_all, reject_death_all, accept_move_all, reject_move_all
!f2py intent(out) accept_sigma_y_all, reject_sigma_y_all

 
 !    call OMP_SET_NUM_THREADS(nbeta)     ! UNCOMMENT IF PARALLEL TEMPERING REQUIRED

  call init_random_seed()          ! Ensure random number generation starts from new point each time program is run
                                  ! Random seed only needs to be called once in a program.  
				    ! It is used to generate the first random number. 
				   ! Any random numbers generated after the first will use the previous random number as its seed.

! SET-UP INITIAL STARTING VALUES

accept(:)=0
accept_birth=0
accept_death=0
accept_move=0
reject(:)=0
reject_birth=0
reject_death=0
reject_move=0
accept_sigma_y=0
reject_sigma_y=0
accept_swap=0
reject_swap=0
it_sub=1

accept_all(:,:)=0
accept_birth_all(:)=0
accept_death_all(:)=0
accept_move_all(:)=0
reject_all(:,:)=0
reject_birth_all(:)=0
reject_death_all(:)=0
reject_move_all(:)=0
accept_sigma_y_all(:)=0
reject_sigma_y_all(:)=0

acc_h_batch(:)=0.
rej_h_batch(:)=0.
accept_batch(:)=0
reject_batch(:) = 0
acc_y_batch(:)=0
rej_y_batch(:)=0

 do jj=1,ydim2   
        y_error_temp(R_indices(:,jj)) = sigma_model(jj,1)    ! Provided dim2 isn't too big then should be fine
 enddo  
   			   
 sigma_y_temp=sqrt(y_error_temp**2 + sigma_measure**2)      ! sigma_y is the combination of model and measurement error
 n0T(:)=sum((n0(:,1)/sigma_y_temp)**2)

 do ibeta = 1,nbeta
    sigma_y(:,ibeta) = sigma_y_temp
 enddo
 sigma_model_ap=sigma_model(:,1)
 detval(:) = sum(alog(sigma_y(:,1)))


! MCMC loop
!###############################################################
do it=1,(nIt+burn_in)
   
   !call random_number(u) 

   if (rjmcmc .EQ. 1) then    
       ! If doing reversible jump:
       !remain_it = FLOOR(5*u) + 1    ! Choose random number between 1 and 5 to chose what to update
       remain_it= modulo(it,5)+1
   else
       ! If doing fixed dimension MCMC:
       !remain_it = FLOOR(2*u) + 1    ! Choose random number between 1 and 2 - no reversible jump.
       remain_it = modulo(it,2)+1
   endif


!$OMP PARALLEL DO DEFAULT(SHARED) private(ibeta, betaib,kib,xib,pdf_param1ib,pdf_param2ib, plonib, platib), &
!$OMP& private(regions_vib,h_aggib,n0ib,n0Tib,kIC,xib1,n0ib1,n0Tib1,acceptib1,rejectib1,regions_vib1), &
!$OMP& private(u, kib1, h_aggib1, plonib1, platib1, acceptxib1, rejectxib1), &
!$OMP& private(sigma_yib, sigma_modelib, sigma_yib1, sigma_modelib1, accept_yib1, reject_yib1), &
!$OMP& private(pdf_param1ib1, pdf_param2ib1, detvalib, detvalib1),&
!$OMP& private(stepsize_ib1, acc_bxib1, rej_bxib1),&
!$OMP& private(stepsize_p1_ib1, acc_prob_p1_ib1,rej_prob_p1_ib1),&
!$OMP& private(stepsize_p2_ib1, stepsize_sig_ib1, acc_byib1, rej_byib1),&
!$OMP& shared(x,n0,n0T, k, pdf_param1, pdf_param2, h_agg, plon,plat, regions_v)
   do ibeta=1,nbeta


     if (para_temp .EQ. 1 .or. ibeta .EQ. 1) then 

       ! The following ib variables are necessary for PT
       betaib = beta(ibeta)
       kib = k(ibeta)
       xib  = x(:,ibeta)
       pdf_param1ib = pdf_param1(:,ibeta)
       pdf_param2ib = pdf_param2(:,ibeta)

       plonib = plon(:,ibeta)
       platib = plat(:,ibeta)
       regions_vib = regions_v(:,ibeta)
       h_aggib = h_agg(:,:,ibeta)
       n0ib = n0(:,ibeta)
       n0Tib = n0T(ibeta)

       sigma_yib = sigma_y(:,ibeta)
       sigma_modelib = sigma_model(:,ibeta)
       detvalib = detval(ibeta)

       kIC = kib+nIC
       
       
       if (remain_it .EQ. 1) then              ! X UPDATE AND EMISSIONS HYPERPARAMETER UPDATE

            call x_hparam_update(betaib, kib, xib, pdf_param1ib,pdf_param2ib, &
                 pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf, &
                 pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf, &
                 acc_h_batch, rej_h_batch, x_pdf_all, it, burn_in, nIC, kICmax, nIC1, &
                 pdf_param1ib1, pdf_param2ib1, stepsize_p1_ib1, stepsize_p2_ib1, acc_prob_p1_ib1, rej_prob_p1_ib1) 


            call x_update(betaib,kib, xib, pdf_param1ib1,pdf_param2ib1, &
                         h_aggib,n0ib,n0Tib,sigma_yib, stepsize, &
                         accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
                         xib1, n0ib1, n0Tib1, acceptxib1, rejectxib1, stepsize_ib1, acc_bxib1, rej_bxib1)


            x(:,ibeta) = xib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 
            pdf_param1(:,ibeta) = pdf_param1ib1
            pdf_param2(:,ibeta) = pdf_param2ib1
            !stepsize=stepsize_ib1
            !stepsize_pdf_p1=stepsize_p1_ib1
            !stepsize_pdf_p2=stepsize_p2_ib1
            !acc_h_batch=acc_prob_p1_ib1
            !rej_h_batch=rej_prob_p1_ib1
            !accept_batch=acc_bxib1
            !reject_batch=rej_bxib1

            accept_all(:,ibeta) = accept_all(:,ibeta) + acceptxib1
            reject_all(:,ibeta) = reject_all(:,ibeta) + rejectxib1

            if (betaib .EQ. 1.) then 
               accept(:) = accept(:) + acceptxib1
               reject(:) = reject(:) + rejectxib1
               stepsize=stepsize_ib1
               accept_batch=acc_bxib1
               reject_batch=rej_bxib1
               stepsize_pdf_p1=stepsize_p1_ib1
               stepsize_pdf_p2=stepsize_p2_ib1
               acc_h_batch=acc_prob_p1_ib1
               rej_h_batch=rej_prob_p1_ib1
            endif
            
           

       elseif (remain_it .EQ. 3) then       ! BIRTH
               
               call birth(betaib,kib, xib, h_aggib,y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lons,lats, & 
                          h_v, pdf_param1ib, pdf_param2ib, x_pdf_all(nIC1), &
                          lonmin, lonmax, &
                          latmin,latmax, sigma_bd,it,burn_in,nIC,kICmax,kmax, &
                          nmeasure,Ngrid, &
                          kib1, xib1, h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1,&
                          pdf_param1ib1, pdf_param2ib1)

                k(ibeta) = kib1
                x(:,ibeta) = xib1
                plon(:,ibeta) = plonib1
                plat(:,ibeta) = platib1
                regions_v(:,ibeta) = regions_vib1
                h_agg(:,:,ibeta) = h_aggib1
                n0(:,ibeta) = n0ib1
                n0T(ibeta) = n0Tib1
                pdf_param1(:,ibeta) = pdf_param1ib1
                pdf_param2(:,ibeta) = pdf_param2ib1


                accept_birth_all(ibeta)= accept_birth_all(ibeta) + acceptib1
                reject_birth_all(ibeta)= reject_birth_all(ibeta) + rejectib1

               if (betaib .EQ. 1.) then 
                   accept_birth= accept_birth + acceptib1
                   reject_birth= reject_birth + rejectib1
               endif

           elseif (remain_it .EQ. 4) then    ! DEATH

               call death(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lons,lats, & 
                          h_v, pdf_param1ib, pdf_param2ib, x_pdf_all(nIC1), sigma_bd, &
                          it, burn_in,nIC, kICmax, kmin, kmax, nmeasure, &
                          Ngrid, &
                          kib1, xib1, h_aggib1,n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1, &
                          pdf_param1ib1, pdf_param2ib1)
                     
               k(ibeta) = kib1
               x(:,ibeta) = xib1
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               pdf_param1(:,ibeta) = pdf_param1ib1
               pdf_param2(:,ibeta) = pdf_param2ib1

               accept_death_all(ibeta)= accept_death_all(ibeta) + acceptib1
               reject_death_all(ibeta)= reject_death_all(ibeta) + rejectib1
      
               if (betaib .EQ. 1.) then 
                   accept_death= accept_death + acceptib1
                   reject_death= reject_death + rejectib1
               endif

           elseif (remain_it .EQ. 5) then    ! MOVE
                
               
               call move(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lons,lats, & 
                         h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, it, &
                         burn_in, nIC, kICmax, kIC, kmax, nmeasure, Ngrid, &
                         h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)
               
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               
               accept_move_all(ibeta)= accept_move_all(ibeta) + acceptib1
               reject_move_all(ibeta)= reject_move_all(ibeta) + rejectib1

               if (betaib .EQ. 1.) then 
                  accept_move= accept_move + acceptib1
                  reject_move= reject_move + rejectib1
               endif

          elseif (remain_it .EQ. 2) then  ! SIGMA_Y UPDATE
             
              call sigma_y_update(betaib, sigma_modelib, sigma_model_ap, sigma_measure, sigma_yib, detvalib, &
                 sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, R_indices, &
                 n0ib,n0Tib, acc_y_batch, rej_y_batch, it, burn_in, nmeasure, ydim1, ydim2, &
                 n0Tib1, accept_yib1, reject_yib1, sigma_yib1, sigma_modelib1, detvalib1, &
                 stepsize_sig_ib1, acc_byib1, rej_byib1) 

              sigma_y(:,ibeta) = sigma_yib1
              sigma_model(:,ibeta) = sigma_modelib1
              n0T(ibeta) = n0Tib1
              detval(ibeta) = detvalib1 

              accept_sigma_y_all(ibeta) = accept_sigma_y_all(ibeta) + accept_yib1
              reject_sigma_y_all(ibeta) = reject_sigma_y_all(ibeta) + reject_yib1

              if (betaib .EQ. 1.) then 
               accept_sigma_y = accept_sigma_y + accept_yib1
               reject_sigma_y = reject_sigma_y + reject_yib1
               stepsize_sigma_y=stepsize_sig_ib1
               acc_y_batch=acc_byib1
               rej_y_batch=rej_byib1
              endif

           endif     ! remain_it
           
      endif     !para_temp .EQ. 1 .or. ibeta .EQ. 1) 

   enddo    ! beta loop
!$OMP END PARALLEL DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Store xit and swap betas
  

   IF (it .GT. burn_in/2) THEN      ! Begin swaps after half of burn-in time
        remain_swap = modulo(it,2)
        IF (remain_swap .EQ. 1) THEN
          if (para_temp .EQ. 1) then
            call random_number(u1)   
            pair1 = FLOOR(nbeta*u1) + 1
            call random_number(u2)  
            pair2 = FLOOR(nbeta*u2) + 1
            if (pair1 .EQ. pair2) THEN
                if (pair2 .EQ. 1) then
                    pair2=pair1+1
                else 
                    pair2=pair1-1
                endif
            endif

            beta1=beta(pair1)*1.
            beta2=beta(pair2)*1.
            pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+detval(pair2)-detval(pair1))  ! detvals should be inverse determinants so signs this way round
            call random_number(randomu)
            if (alog(randomu) .LE. pT_chain) then           
                    beta(pair2)=beta1*1.         ! Uncomment for PT
                    beta(pair1)=beta2*1.         ! Uncomment for PT
                    accept_swap=accept_swap+1              
            else
                reject_swap=reject_swap+1
            endif      ! pT_chain if   
           else
                reject_swap=reject_swap+1
           endif   ! para_temp=1  
         ENDIF      ! reamin_swap =0 if
   ENDIF          ! it > burn_in/2
   


   IF (it .GT. burn_in) THEN     
        remain = modulo(it,nsub)          ! nsub typically = 100
        if (remain .EQ. 0) then

           do ib=1,nbeta                             ! uncomment for PT
               if (beta(ib) .EQ. 1.) then            ! uncomment for PT
          !        ib=1                   ! DEFAULT - AGAIN ib=1 comment and uncomment if statement if doing PT
                  ! STORE THE FOLLOWING VARIABLES AT THINNED nsub FREQUENCY
                  x_it(:,it_sub)=x(:,ib)
                  plon_it(:,it_sub)=plon(:,ib)
                  plat_it(:,it_sub)=plat(:,ib)
                  k_it(it_sub)=k(ib)
                  regions_it(:,it_sub)=regions_v(:,ib)
                  sigma_model_it(:,it_sub)=sigma_model(:,ib)
                  sigma_y_it(:,it_sub)=sigma_y(:,ib)
                  n0T_it(it_sub)=n0T(ib)
                  pdf_param1_it(:,it_sub)=pdf_param1(:,ib)
                  pdf_param2_it(:,it_sub)=pdf_param2(:,ib)     
                  y_it(:,it_sub) = matmul(h_agg(:,:,ib), x(:,ib))  
                  it_sub=it_sub+1
               endif                         ! uncomment for PT
            enddo                            ! uncomment for PT
        endif
   ENDIF           ! it >= burn_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo  ! It loop    ! END OF MCMC LOOP

! OUTPUTS
k_out = k_it
x_out=x_it
y_out = y_it
regions_out = regions_it
plon_out = plon_it
plat_out=plat_it
sigma_model_out=sigma_model_it
sigma_y_out=sigma_y_it
n0T_out=n0T_it
pdf_param1_out=pdf_param1_it
pdf_param2_out=pdf_param2_it
tot_acc_sigma_y=stepsize_sigma_y
tot_acc_x=stepsize
tot_acc_p1=stepsize_pdf_p1
tot_acc_p2=stepsize_pdf_p2
END SUBROUTINE hbtdmcmc


SUBROUTINE x_hparam_update(beta, k, x, pdf_param1_all,pdf_param2_all, &
pdf_p1_hparam1_all, pdf_p1_hparam2_all, stepsize_pdf_p1, pdf_param1_pdf, &
pdf_p2_hparam1_all, pdf_p2_hparam2_all, stepsize_pdf_p2, pdf_param2_pdf, &
accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nIC1, &
pdf_param1_out, pdf_param2_out, stepsize_p1_out, stepsize_p2_out, accept_batch_out, reject_batch_out) 


Implicit none
INTEGER k, nIC, kICmax, nIC1, burn_in, it
REAL av_acc, beta
REAL x(kICmax) 
INTEGER x_pdf_all(nIC1)
INTEGER accept_batch(nIC1), reject_batch(nIC1)
REAL pdf_param1_all(kICmax), pdf_param2_all(kICmax)
REAL pdf_p1_hparam1_all(nIC1), pdf_p1_hparam2_all(nIC1) 
REAL  pdf_p2_hparam1_all(nIC1), pdf_p2_hparam2_all(nIC1)
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)
REAL stepsize_p1_out(nIC1), stepsize_p2_out(nIC1)
INTEGER accept_batch_out(nIC1), reject_batch_out(nIC1)
INTEGER elem(5)
REAL u(5)
REAL pdf_param1        
REAL pdf_param2
REAL accep_prob(nIC1)
!INTEGER accept(nIC1), reject(nIC1)
REAL stepsize_pdf_p2(nIC1), stepsize_pdf_p1(nIC1)
INTEGER xi, x_pdf, xx
REAL pT, randomu, p0,p1
INTEGER pdf_param2_pdf, pdf_param1_pdf
REAL pdf_p2_hparam1, pdf_p2_hparam2, pdf_p1_hparam1, pdf_p1_hparam2
REAL dpdf_param2, pdf_param2_new, dpdf_param1, pdf_param1_new
REAL p0_pdf_param2, p1_pdf_param2, p0_pdf_param1, p1_pdf_param1
REAL stepsize_pdf_p10, stepsize_pdf_p20

! Propose new values of pdf_param1 and pdf_param2

   call random_number(u)   
    
  elem = FLOOR(k*u)+1+nIC

  do xx=1,nIC+5
  
  if (xx .LE. nIC) then
     xi = xx
  else
     xi=elem(xx-nIC)
  endif

  if (xi .LE. nIC) then
     x_pdf = x_pdf_all(xi)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(xi)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(xi)
     pdf_p2_hparam1 = pdf_p2_hparam1_all(xi)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(xi)
     stepsize_pdf_p10 = stepsize_pdf_p1(xi)
     stepsize_pdf_p20 = stepsize_pdf_p2(xi)
  else if (xi .GT. nIC) then
     x_pdf = x_pdf_all(nIC1)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(nIC1)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(nIC1)
     pdf_p2_hparam1 = pdf_p2_hparam1_all(nIC1)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(nIC1) 
     stepsize_pdf_p10 = stepsize_pdf_p1(nIC1)
     stepsize_pdf_p20 = stepsize_pdf_p2(nIC1)
  endif
  dpdf_param1 = random_normal()*stepsize_pdf_p10
  pdf_param1_new = pdf_param1 + dpdf_param1

  dpdf_param2 = random_normal()*stepsize_pdf_p20
  pdf_param2_new = pdf_param2 + dpdf_param2
  
         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior

         call calc_pdf(pdf_param1,pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p0_pdf_param1) 
         call  calc_pdf(pdf_param1_new, pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p1_pdf_param1)

         call calc_pdf(pdf_param2,pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p0_pdf_param2) 
         call  calc_pdf(pdf_param2_new, pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p1_pdf_param2)

      
         call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(x(xi),pdf_param1_new,pdf_param2_new,x_pdf, p1)        
               
         pT = p1+p1_pdf_param1+p1_pdf_param2-p0-p0_pdf_param1-p0_pdf_param2  ! All p1s already in log space
           ! Likelihood values will remain unchanged since x doesn't change

         if (pdf_param2_pdf .eq. 1) then
             if (pdf_param2_new .lt. pdf_p2_hparam1) pT = -1.e20
             if (pdf_param2_new .gt. pdf_p2_hparam2) pT = -1.e20
         endif
         if (pdf_param1_pdf .eq. 1) then
             if (pdf_param1_new .lt. pdf_p1_hparam1) pT = -1.e20
             if (pdf_param1_new .gt. pdf_p1_hparam2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN
             !ACCEPT   
             pdf_param1_all(xi)=pdf_param1_new
             pdf_param2_all(xi)=pdf_param2_new  

             if(beta .eq. 1. .and. it .le. burn_in) then
             !if(it .le. burn_in) then  
                 if (xi .LE. nIC) then
                    accept_batch(xi) = accept_batch(xi) + 1
                 else if (xi .GT. nIC) then
                    accept_batch(nIC1) = accept_batch(nIC1) + 1
                 endif
             endif      ! it le burn_in

         else
             if(beta .eq. 1. .and. it .le. burn_in) then  
                if (xi .LE. nIC) then 
                   reject_batch(xi) = reject_batch(xi) + 1
                else
                   reject_batch(nIC1) = reject_batch(nIC1) + 1  
                endif   ! xi le nIC
             endif    ! it le burn_in
         endif   ! randomu condition

         if(beta .eq. 1.) then
         if(it .le. burn_in .and. modulo(it,500) .eq. 0) then
         !if(it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100*nIC1) then
             if (xi .LE. nIC) then
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) then
                         stepsize_pdf_p1(xi) = exp(alog(stepsize_pdf_p1(xi)) - av_acc)
                         stepsize_pdf_p2(xi) = exp(alog(stepsize_pdf_p2(xi)) - av_acc)
                     endif
                     if(accep_prob(xi) .gt. 0.6) then
                         stepsize_pdf_p1(xi) = exp(alog(stepsize_pdf_p1(xi)) + av_acc)
                         stepsize_pdf_p2(xi) = exp(alog(stepsize_pdf_p2(xi)) + av_acc)
                     endif
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
             else if (xi .GT. nIC) then
                 if(accept_batch(nIC1)+reject_batch(nIC1) .gt. 0) then    
                      accep_prob(nIC1) = real(accept_batch(nIC1))/(accept_batch(nIC1)+reject_batch(nIC1))
                      av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(nIC1)+reject(nIC1)))
                      if(accep_prob(nIC1) .lt. 0.2) then
                         stepsize_pdf_p1(nIC1) = exp(alog(stepsize_pdf_p1(nIC1)) - av_acc)
                         stepsize_pdf_p2(nIC1) = exp(alog(stepsize_pdf_p2(nIC1)) - av_acc)
                      endif
                      if(accep_prob(nIC1) .gt. 0.6) then
                         stepsize_pdf_p1(nIC1) = exp(alog(stepsize_pdf_p1(nIC1)) + av_acc)
                         stepsize_pdf_p2(nIC1) = exp(alog(stepsize_pdf_p2(nIC1)) + av_acc)
                      endif
                      accept_batch(nIC1) = 0
                      reject_batch(nIC1) = 0
                 endif
             endif

          endif    ! modulo(it,500)
          endif    ! beta=1
  enddo


pdf_param1_out=pdf_param1_all
pdf_param2_out=pdf_param2_all
stepsize_p1_out=stepsize_pdf_p1
stepsize_p2_out=stepsize_pdf_p2
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE x_hparam_update

SUBROUTINE x_update(beta,k, x, pdf_param1_all,pdf_param2_all,  &
h_agg,n0,n0T,sigma_y, stepsize,&
accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
x_out, n0_out, n0T_out, accept_out, reject_out, stepsize_out, acc_batch_out, rej_batch_out) 

Implicit none 
INTEGER nmeasure, it, burn_in, k, nIC, kICmax, nIC1
REAL beta, n0T, n0T_out 
REAL av_acc
REAL x(kICmax) 
REAL x_out(kICmax) 
REAL h_agg(nmeasure,kICmax)   
REAL n0(nmeasure) 
REAL n0_out(nmeasure) 
REAL sigma_y(nmeasure)
REAL dy(nmeasure)
REAL n1(nmeasure) 
INTEGER x_pdf_all(nIC1)
REAL pdf_param1_all(kICmax), pdf_param2_all(kICmax)
INTEGER elem(5)
REAL u(5)
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(nIC1)
REAL accep_prob(nIC1)
INTEGER accept(nIC1), reject(nIC1)
INTEGER accept_batch(nIC1)
INTEGER reject_batch(nIC1)
INTEGER accept_out(nIC1), reject_out(nIC1)
REAL stepsize_out(nIC1)
INTEGER xi, x_pdf, xx
REAL dx, n1T, pT, randomu, p0,p1
REAL x_new(kICmax)
REAL stepsize0
INTEGER acc_batch_out(nIC1)
INTEGER rej_batch_out(nIC1)

accept=0
reject=0

! Propose new x values
! Currently this is done for all fixed x values and 5 elements of the variable regions
   ! CHANGE OF EMISSIONS VALUES
  call random_number(u)   
    
  elem = FLOOR(k*u)+1+nIC

  do xx=1,nIC+5
  
  if (xx .LE. nIC) then
     xi = xx
  else
     xi=elem(xx-nIC)
  endif

  if (xi .LE. nIC) then
     x_pdf = x_pdf_all(xi)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     stepsize0 = stepsize(xi)

  else if (xi .GT. nIC) then
     x_pdf = x_pdf_all(nIC1)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     stepsize0 = stepsize(nIC1)
    
  endif
  dx = random_normal()*stepsize0

  x_new = x
  x_new(xi) =x(xi)+dx
  p0=0.
  p1=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
         dy=h_agg(:,xi)*dx 
         n1=n0+dy

         n1T=sum((n1/sigma_y)**2)
          
         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior
                  
         call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(x_new(xi),pdf_param1,pdf_param2,x_pdf, p1)        
       
         pT = p1-p0 -0.5*(n1T - n0T)*beta ! All p1s already in log space

         if (x_pdf .eq. 1) then
             if (x_new(xi) .lt. pdf_param1) pT = -1.e20
             if (x_new(xi) .gt. pdf_param2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN

             !ACCEPT
             x(xi)=x(xi)+dx
             n0=n1
             n0T=n1T
             if (xi .LE. nIC) then
                !if (beta .EQ. 1. .and. it .GT. burn_in) accept(xi) = accept(xi) + 1
                if (it .GT. burn_in) accept(xi) = accept(xi) + 1
             else if (xi .GT. nIC) then
                !if (beta .EQ. 1. .and. it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
                if (it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
             endif
             
             !if (beta .EQ. 1. .and. it .GT. burn_in) accept = accept + 1
                                
         else
             !REJECT
             !if (beta .EQ. 1. .and. it .GT. burn_in) then 
             if (it .GT. burn_in) then 
                 if (xi .LE. nIC) then 
                     reject(xi) = reject(xi) + 1
                 else
                     reject(nIC1) = reject(nIC1) + 1  
                 endif
             endif
         endif   ! randomu condition

         if(beta .eq. 1. .and. it .le. burn_in) then
             if (alog(randomu) .LE. pT) THEN

                 if (xi .LE. nIC) then
                    !if (beta .EQ. 1. .and. it .GT. burn_in) accept_batch(xi) = accept_batch(xi) + 1
                    accept_batch(xi) = accept_batch(xi) + 1
                 else if (xi .GT. nIC) then
                   ! if (beta .EQ. 1. .and. it .GT. burn_in) accept_batch(nIC1) = accept_batch(nIC1) + 1
                    accept_batch(nIC1) = accept_batch(nIC1) + 1
                 endif
             else
                 if (xi .LE. nIC) then 
                     reject_batch(xi) = reject_batch(xi) + 1
                 else
                     reject_batch(nIC1) = reject_batch(nIC1) + 1  
                 endif

             endif

         endif

         if(beta .eq. 1 .and. it .le. burn_in .and. modulo(it,500) .eq. 0) then
         !if(beta .eq. 1 .and. it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100*nIC1) then
             !if (xi .EQ. 1) then
             if (xi .LE. nIC) then
                 !write(*,*) accept_batch
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     !av_acc = min(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) stepsize(xi) = exp(alog(stepsize(xi)) - av_acc)
                     if(accep_prob(xi) .gt. 0.6) stepsize(xi) = exp(alog(stepsize(xi)) + av_acc)
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
             else if (xi .GT. nIC) then
                 if(accept_batch(nIC1)+reject_batch(nIC1) .gt. 0) then    
                      accep_prob(nIC1) = real(accept_batch(nIC1))/(accept_batch(nIC1)+reject_batch(nIC1))
                      !av_acc = min(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(nIC1)+reject(nIC1)))
                      av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                      if(accep_prob(nIC1) .lt. 0.2) stepsize(nIC1) = exp(alog(stepsize(nIC1)) - av_acc)
                      if(accep_prob(nIC1) .gt. 0.6) stepsize(nIC1) = exp(alog(stepsize(nIC1)) + av_acc)
                      accept_batch(nIC1) = 0
                      reject_batch(nIC1) = 0
                 endif
             endif

         endif

  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_out=x
n0_out=n0
n0T_out=n0T
accept_out=accept
reject_out=reject
stepsize_out=stepsize
acc_batch_out = accept_batch
rej_batch_out = reject_batch
END SUBROUTINE x_update



SUBROUTINE birth(beta,k, x, h_agg,y,n0,n0T,sigma_y, plon, plat, regions_v, lons,lats, & 
h_v,pdf_param1, pdf_param2, x_pdf,  &
lonmin, lonmax, latmin,latmax, sigma_bd, &
it, burn_in, nIC, kICmax, kmax, nmeasure, Ngrid, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out, &
pdf_param1_out,pdf_param2_out)

IMPLICIT NONE
! Dimensions
INTEGER nmeasure, Ngrid, k,kmax, nIC, kICmax
REAL lonmin,lonmax,latmin,latmax,sigma_bd
! Single Variables
INTEGER accept_birth, reject_birth, it, burn_in, x_pdf 
REAL beta, n0T
! Input arrays
REAL x(kICmax) 
REAL pdf_param1(kICmax) 
REAL pdf_param2(kICmax)
REAL h_agg(nmeasure,kICmax)  
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL sigma_y(nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lons(Ngrid)
REAL lats(Ngrid)
REAL h_v(nmeasure, Ngrid)
! Outputs
INTEGER k_out
REAL x_out(kICmax) 
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)  
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rib, k1, errstat,jj, kIC      
REAL u, u2, plon_new, plat_new, c_new, x_new, n1Tb, pT_birth, randomu
REAL mu, sigma, pdf_param1_new, pdf_param2_new
INTEGER regions_v1b(Ngrid)       
REAL n1b(nmeasure)     
INTEGER ilon,ilat, reject_stat,zi         
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1b, plat1b, x1b
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2
REAL,PARAMETER     :: pi = 3.14159265 

accept_birth=0
reject_birth=0

! Propose new voronoi cell and increase dimension size of k by 1
k1=k+1
kIC=k1+nIC

if (k1 .LT. kmax) THEN            
  
  allocate(plon1b(k1))     
  allocate(plat1b(k1))     
  allocate(h_agg2(nmeasure, kIC))
  allocate(x1b(kIC))               

   ! 1. Select new nucleus location - different from current locations and on hte underlying grid
  !call random_number(u)
  !plon_new = lonmin+(lonmax-lonmin)*u
  !call random_number(u2)
  !plat_new = latmin+(latmax-latmin)*u2
 
  call random_number(u) 
  call random_number(u2)
  ilon=FLOOR(Ngrid*u)+1   ! Choose random number between 1 and Ngrid
  ilat=FLOOR(Ngrid*u2)+1   ! Choose random number between 1 and Ngrid
        
  plon_new = lons(ilon)
  plat_new = lats(ilat)
  reject_stat=0
  do zi=1,k
     if (plon(zi) .EQ. plon_new .and. plat(zi) .eq. plat_new) then
        reject_stat=1  
     endif
  enddo

  if (reject_stat .ne. 1) then

  plon1b(1:k) = plon(1:k)
  plon1b(k1) = plon_new
  plat1b(1:k) = plat(1:k)
  plat1b(k1) = plat_new
                                  
  ! 2. Recalculate voronoi cells

  call closest_grid(regions_v1b,lons,lats, plon1b(1:k1),plat1b(1:k1), k1, Ngrid)

   h_agg2(:,:)=0.

   if (nIC .GT. 0) then 
       h_agg2(:,1:nIC)=h_agg(:,1:nIC)
   endif
       
    do ri=1,k1
       do jj=1,Ngrid
         if (regions_v1b(jj) .EQ. ri) then
            h_agg2(:,ri+nIC) = h_agg2(:,ri+nIC)+h_v(:,jj)
         endif
      enddo
   enddo
 
 !#######################################################################
                    
  call closest_grid_sng(rib,plon_new,plat_new, plon(1:k),plat(1:k), k) 
                                  
  c_new = x(rib+nIC)
  x_new = random_normal()*sigma_bd+c_new

  pdf_param1_new = pdf_param1(rib+nIC)
  pdf_param2_new = pdf_param2(rib+nIC)
  
  if (x_pdf .EQ. 1) THEN  ! 1=UNIFORM  
                
      if (x_new .GT. pdf_param1_new .and. x_new .LT. pdf_param2_new) THEN
                       
          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new

          n1b=matmul(h_agg2,x1b)-y                    
          n1Tb=sum((n1b/sigma_y)**2)

          pT_birth = alog(sqrt(2*pi)*sigma_bd/(pdf_param2_new-pdf_param1_new)) &
              + ((x_new-c_new)**2)/2./sigma_bd**2 - 0.5*(n1Tb - n0T)*beta
       

         call random_number(randomu)             
         if (alog(randomu) .LE. pT_birth) THEN
           !ACCEPT
           k=k1
           x(:)=0.
           x(1:kIC)=x1b(1:kIC)
           pdf_param1(kIC)=pdf_param1_new
           pdf_param2(kIC)=pdf_param2_new
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1   
           if (it .GT. burn_in) accept_birth=accept_birth+1    
         else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1 
           if (it .GT. burn_in) reject_birth=reject_birth+1  
         endif

      else
          !REJECT if x_new is negative
          !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1
          if (it .GT. burn_in) reject_birth=reject_birth+1
      endif   

  else if (x_pdf .GE. 2) THEN     ! 2=GAUSSIAN 3=LOGNORMAL 

          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new
          
          n1b=matmul(h_agg2,x1b)-y                                         
          n1Tb=sum((n1b/sigma_y)**2)

          if (x_pdf .EQ. 2) THEN   !2=GAUSSIAN

              pT_birth = alog(sigma_bd/pdf_param2_new) + ((x_new-c_new)**2)/2./sigma_bd**2 - &
                         ((x_new-pdf_param1_new)**2)/2./pdf_param2_new**2 &
                         - 0.5*(n1Tb - n0T)*beta  

          else if (x_pdf .EQ. 3) THEN
 
              mu = alog(pdf_param1_new) - 0.5*alog(1. + pdf_param2_new**2/pdf_param1_new**2)
              sigma=sqrt(alog((pdf_param2_new/pdf_param1_new)**2 + 1.))

               pT_birth = alog(sigma_bd/sigma/x_new) + ((x_new-c_new)**2)/2./sigma_bd**2 - &
                         ((alog(x_new)-mu)**2)/2./sigma**2 &
                          - 0.5*(n1Tb - n0T)*beta             

         else
             stop
         endif
       
         call random_number(randomu)             
         if (alog(randomu) .LE. pT_birth) THEN
           !ACCEPT
           k=k1
           x(:)=0.
           x(1:kIC)=x1b(1:kIC)
           pdf_param1(kIC)=pdf_param1_new
           pdf_param2(kIC)=pdf_param2_new
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1    
           if (it .GT. burn_in) accept_birth=accept_birth+1    
         else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1  
           if (it .GT. burn_in) reject_birth=reject_birth+1  
         endif

  endif       ! x_pdf
 else
    !REJECT if plon_new and plat_new are the same as any other point
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1 
    if (it .GT. burn_in) reject_birth=reject_birth+1  
 endif
        
else
    !REJECT if k1 > kmax
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1 
    if (it .GT. burn_in) reject_birth=reject_birth+1  
endif


!! Deallocate arrays in each loop
if(allocated(plon1b))  deallocate(plon1b,stat=errstat)
  if (errstat /= 0) stop
if(allocated(plat1b))  deallocate(plat1b,stat=errstat)
  if (errstat /= 0) stop
if(allocated(h_agg2))  deallocate(h_agg2,stat=errstat)
  if (errstat /= 0) stop
if(allocated(x1b))  deallocate(x1b,stat=errstat)
  if (errstat /= 0) stop

k_out=k
x_out=x
pdf_param1_out=pdf_param1
pdf_param2_out=pdf_param2
h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_birth
reject_out=reject_birth


END SUBROUTINE birth


SUBROUTINE death(beta,k, x, h_agg,y,n0,n0T,sigma_y, plon, plat, regions_v, lons,lats, & 
h_v, pdf_param1, pdf_param2, x_pdf, &
sigma_bd, &
it, burn_in, nIC, kICmax, kmin, kmax, nmeasure, Ngrid, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, &
plon_out, plat_out, accept_out, reject_out, pdf_param1_out, pdf_param2_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax,nmeasure,Ngrid, accept_death, reject_death, it, burn_in, kmin, k, nIC, kICmax
REAL sigma_bd, beta, n0T 
! Input arrays
REAL x(kICmax) 
REAL pdf_param1(kICmax)
REAL pdf_param2(kICmax)
REAL h_agg(nmeasure,kICmax) 
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL sigma_y(nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lons(Ngrid)
REAL lats(Ngrid)
REAL h_v(nmeasure, Ngrid)
INTEGER x_pdf
! Outputs
INTEGER k_out 
REAL n0T_out
REAL x_out(kICmax)
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rid, k1d, jj, ci_rm, errstat, kIC     
REAL u, plon_rm, plat_rm, x_cell, x_rm, n1Td, pT_death, randomu
REAL mu,sigma, pdf_param1_rm, pdf_param2_rm
INTEGER regions_v1d(Ngrid)      
REAL n1d(nmeasure)             
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1d, plat1d, x1d, pdf_param1d, pdf_param2d 
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2d

REAL,PARAMETER     :: pi = 3.14159265 

accept_death=0
reject_death=0

!DEATH
k1d=k-1

kIC=k1d+nIC

if (k1d .GE. kmin) THEN            

  allocate(plon1d(k1d))     
  allocate(plat1d(k1d))     
  allocate(h_agg2d(nmeasure, kIC))
  allocate(x1d(kIC)) 
  allocate(pdf_param1d(kIC)) 
  allocate(pdf_param2d(kIC))  
           
  ! 1. Select new cell location - needs to be different to current locations
   
  call random_number(u)   
  ci_rm = FLOOR((k1d+1)*u) + 1
  
  plon_rm = plon(ci_rm)
  plat_rm = plat(ci_rm)
  x_rm = x(ci_rm+nIC)

  pdf_param1_rm = pdf_param1(ci_rm+nIC)  
  pdf_param2_rm = pdf_param2(ci_rm+nIC) 

  IF (nIC .GT. 0) THEN
     x1d(1:nIC) = x(1:nIC)
     pdf_param1d(1:nIC) = pdf_param1(1:nIC)
     pdf_param2d(1:nIC) = pdf_param2(1:nIC)
  ENDIF

  IF (ci_rm .EQ. 1) THEN
      plon1d(1:k1d) = plon(2:(k1d+1))
      plat1d(1:k1d) = plat(2:(k1d+1))
      x1d(nIC+1:kIC) = x(nIC+2:(kIC+1))
      pdf_param1d(nIC+1:kIC) = pdf_param1(nIC+2:(kIC+1))  
      pdf_param2d(nIC+1:kIC) = pdf_param2(nIC+2:(kIC+1))  

  ELSEIF (ci_rm .EQ. (k1d+1)) THEN
      plon1d(1:k1d) = plon(1:k1d)
      plat1d(1:k1d) = plat(1:k1d)
      x1d(nIC+1:kIC) = x(nIC+1:kIC)
      pdf_param1d(nIC+1:kIC) = pdf_param1(nIC+1:kIC)  
      pdf_param2d(nIC+1:kIC) = pdf_param2(nIC+1:kIC)  

  ELSE   
      plon1d(1:(ci_rm-1)) = plon(1:(ci_rm-1))                   
      plon1d(ci_rm:k1d) = plon((ci_rm+1):(k1d+1))
      plat1d(1:(ci_rm-1)) = plat(1:(ci_rm-1))                   
      plat1d(ci_rm:k1d) = plat((ci_rm+1):(k1d+1)) 

      x1d(nIC+1:(ci_rm+nIC-1)) = x(nIC+1:(ci_rm+nIC-1)) 
      x1d(ci_rm+nIC:kIC) = x((ci_rm+nIC+1):(kIC+1))
      pdf_param1d(nIC+1:(ci_rm+nIC-1)) = pdf_param1(nIC+1:(ci_rm+nIC-1)) 
      pdf_param1d(ci_rm+nIC:kIC) = pdf_param1((ci_rm+nIC+1):(kIC+1))
      pdf_param2d(nIC+1:(ci_rm+nIC-1)) = pdf_param2(nIC+1:(ci_rm+nIC-1)) 
      pdf_param2d(ci_rm+nIC:kIC) = pdf_param2((ci_rm+nIC+1):(kIC+1))

  ENDIF
                                            
  ! 2. Recalculate voronoi cells

  call closest_grid(regions_v1d,lons,lats, plon1d,plat1d, k1d, Ngrid)

   h_agg2d(:,:)=0.
  
   if (nIC .GT. 0) then
       h_agg2d(:,1:nIC)=h_agg(:,1:nIC)
   endif

   do ri=1,k1d
      do jj=1,Ngrid
         if (regions_v1d(jj) .EQ. ri) then
            h_agg2d(:,ri+nIC) = h_agg2d(:,ri+nIC)+h_v(:,jj)
         endif
      enddo
   enddo

 !#######################################################################
                    
  call closest_grid_sng(rid,plon_rm,plat_rm, plon1d,plat1d, k1d) 
                                  
  x_cell = x1d(rid+nIC)
                    
  n1d=matmul(h_agg2d,x1d)-y
  n1Td=sum((n1d/sigma_y)**2)                           
  ! ACCEPTANCE PROBABILITY  

  IF (x_pdf .EQ. 1) THEN  ! 1 = UNIFORM
     pT_death = alog((pdf_param2_rm-pdf_param1_rm)/sqrt(2.*pi)/sigma_bd) &
      - ((x_cell-x_rm)**2)/2./sigma_bd**2 -0.5*(n1Td - n0T)*beta

  ELSE IF (x_pdf .EQ. 2) THEN   ! 2=GAUSSIAN

      pT_death = alog(pdf_param2_rm/sigma_bd) - ((x_cell-x_rm)**2)/2./sigma_bd**2 + &
                         ((x_rm-pdf_param1_rm)**2)/2./pdf_param2_rm**2 &
                          - 0.5*(n1Td - n0T)*beta


  ELSE IF (x_pdf .EQ. 3) THEN   ! 3=LOGNORMAL

      mu = alog(pdf_param1_rm) - 0.5*alog(1. + pdf_param2_rm**2/pdf_param1_rm**2)
      sigma=sqrt(alog((pdf_param2_rm/pdf_param1_rm)**2 + 1.))
  
      pT_death = alog(sigma*x_rm/sigma_bd) - ((x_cell-x_rm)**2)/2./sigma_bd**2 + &
                         ((alog(x_rm)-mu)**2)/2./sigma**2 &
                          - 0.5*(n1Td - n0T)*beta

   
  ELSE
     ! write(*,*) 'There is an error with x_pdf in death loop - exiting!'
      stop
  ENDIF 

  !write(*,*) pT_death                             
  call random_number(randomu)             
  if (alog(randomu) .LE. pT_death) THEN
          !ACCEPT
           k=k1d
           x(:)=0.
           x(1:kIC)=x1d
           pdf_param1(:)=0.
           pdf_param1(1:kIC)=pdf_param1d
           pdf_param2(:)=0.
           pdf_param2(1:kIC)=pdf_param2d
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2d
           n0=n1d
           n0T=n1Td
           regions_v(:)=regions_v1d
           plon(:)=0.
           plat(:)=0.
           plon(1:k1d)=plon1d
           plat(1:k1d)=plat1d
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_death=accept_death+1   
           if (it .GT. burn_in) accept_death=accept_death+1   

       else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1  
           if (it .GT. burn_in) reject_death=reject_death+1  
       endif
        
else
    !REJECT if k1d < kmin
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1 
    if (it .GT. burn_in) reject_death=reject_death+1 
endif

!! Deallocate arrays in each loop
if(allocated(plon1d))  deallocate(plon1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(plat1d))  deallocate(plat1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(h_agg2d))  deallocate(h_agg2d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(x1d))  deallocate(x1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(pdf_param1d))  deallocate(pdf_param1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(pdf_param2d))  deallocate(pdf_param2d,stat=errstat)
  if (errstat /= 0) stop


k_out=k
x_out=x
pdf_param1_out=pdf_param1
pdf_param2_out=pdf_param2
h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_death
reject_out=reject_death


END SUBROUTINE death


SUBROUTINE move(beta,k, x, h_agg, y,n0,n0T,sigma_y, plon, plat, regions_v, lons,lats, & 
h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, it, &
burn_in, nIC, kICmax, kIC, kmax, nmeasure, Ngrid, &
h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax, nmeasure, Ngrid, accept_move, reject_move, it, burn_in,k, nIC, kIC, kICmax
! Single Variables
REAL lonmin,lonmax,latmin,latmax, sigma_clon, sigma_clat, beta 
! Input arrays
REAL x(kICmax) 
REAL h_agg(nmeasure,kICmax)  
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL n0T
REAL sigma_y(nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lons(Ngrid)
REAL lats(Ngrid)
REAL h_v(nmeasure, Ngrid)
REAL plon1m(k)
REAL plat1m(k)
REAL x1m(kIC)
REAL h_agg2m(nmeasure,kIC)
! Outputs
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  ri, k1, jj, ci_mv    
REAL u, n1Tm, pT_move, randomu
INTEGER regions_v1m(Ngrid)      
REAL n1m(nmeasure) 
INTEGER reject_stat, zi
! Allocatable arrays
! None
REAL,PARAMETER     :: pi = 3.14159265 

accept_move=0
reject_move=0

   !MOVE
   k1=k
           
   ! 1. Select new cell location - needs to be different to current locations

   x1m=x(1:kIC)               
   plon1m = plon(1:k1)
   plat1m = plat(1:k1)
   call random_number(u)   
   ci_mv = FLOOR((k1)*u) + 1
                
   ! 2. Move voronoi cell to new random location
   ! Location based on gaussian PDF about current position
    
   !plon1m(ci_mv) = random_normal()*sigma_clon+plon1m(ci_mv)
   !plat1m(ci_mv) = random_normal()*sigma_clat+plat1m(ci_mv)       

   ! Move onto centre of underlying grid 
   plon1m(ci_mv) = FLOOR(random_normal()*sigma_clon)*(lons(2)-lons(1))+plon1m(ci_mv)
   plat1m(ci_mv) = FLOOR(random_normal()*sigma_clat)*(lats(2)-lats(1))+plat1m(ci_mv)
            
   reject_stat=0
   do zi=1,k
     if (plon(zi) .EQ. plon1m(ci_mv) .and. plat(zi) .eq. plat1m(ci_mv)) then
        reject_stat=1  
     endif
   enddo    

   if (reject_stat .ne. 1) then      
         
   ! Need to reject if outside of lon/lat range.
   IF (plon1m(ci_mv) .GT. lonmax .OR. plon1m(ci_mv) .LT. lonmin) THEN
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1
       if (it .GT. burn_in) reject_move=reject_move+1
                         
   ELSEIF (plat1m(ci_mv) .GT. latmax .OR. plat1m(ci_mv) .LT. latmin) THEN           
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1 
       if (it .GT. burn_in) reject_move=reject_move+1

   ELSE    
                                              
      ! 2. Recalculate voronoi cells

      call closest_grid(regions_v1m,lons,lats, plon1m,plat1m, k1, Ngrid)

      h_agg2m(:,:)=0.
      
      if (nIC .GT. 0) then 
          h_agg2m(:,1:nIC)=h_agg(:,1:nIC)
      endif

      do ri=1,k1
         do jj=1,Ngrid
            if (regions_v1m(jj) .EQ. ri) then
                h_agg2m(:,ri+nIC) = h_agg2m(:,ri+nIC)+h_v(:,jj)
            endif
         enddo
      enddo

 !#######################################################################
     
     n1m=matmul(h_agg2m,x1m)-y
     n1Tm=sum((n1m/sigma_y)**2)
                                     
     ! ACCEPTANCE PROBABILITY   
     pT_move = (n1Tm - n0T)*(-0.5)*beta  
                              
     call random_number(randomu)             
     if (alog(randomu) .LE. pT_move) THEN
           !ACCEPT
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2m
           n0=n1m
           n0T=n1Tm
           regions_v=regions_v1m
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1m
           plat(1:k1)=plat1m
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_move=accept_move+1   
           if (it .GT. burn_in) accept_move=accept_move+1   
      else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1  
           if (it .GT. burn_in) reject_move=reject_move+1
      endif
        
   
  ENDIF    ! plon & lat within range
   
else   ! lon_new, lat_new on same location as another point
      !REJECT
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1  
       if (it .GT. burn_in) reject_move=reject_move+1
endif     

h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_move
reject_out=reject_move

END SUBROUTINE move

SUBROUTINE sigma_y_update(beta, sigma_model_current, sigma_model_ap, sigma_measure, sigma_y_current, &
detval_current,sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, R_indices, &
n0,n0T, accept_batch, reject_batch, it, burn_in, nmeasure, dim1, dim2, &
n0T_out, accept_out, reject_out, sigma_y_out, sigma_model_out, detval_out, &
stepsize_sig_out, accept_batch_out, reject_batch_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, dim1, dim2
! Single Variables
REAL beta 
! Input arrays
REAL n0(nmeasure) 
REAL n0T
REAL detval_current
REAL sigma_measure(nmeasure)
REAL sigma_model_current(dim2)
REAL sigma_model_ap(dim2)
REAl sigma_y_current(nmeasure)
INTEGER R_indices(dim1,dim2)
REAL stepsize_sigma_y(dim2)
REAL sigma_model_hparam1(dim2)
REAL sigma_model_hparam2(dim2)
INTEGER sigma_model_pdf
! Outputs
REAL n0T_out, detval_out
REAL sigma_y_out(nmeasure)
REAL sigma_model_out(dim2)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  yi, jj,ii
REAL randomu, dsigma_y, sigma_model_new, u, av_acc
REAL p0_sigma_y, p1_sigma_y, n1T, detval_new, pT
REAL y_error_new(nmeasure)
REAL sigma_y_new(nmeasure)
REAL stepsize_sig_out(dim2)
INTEGER accept_batch(dim2), reject_batch(dim2)
INTEGER accept_batch_out(dim2), reject_batch_out(dim2)
REAL accep_prob

accept=0
reject=0

 !do yi=1,dim2
    call random_number(u)   
    yi = FLOOR(dim2*u)+1
		
    ! Generate new value of sigma_y
    dsigma_y = random_normal()*stepsize_sigma_y(yi)*sigma_model_ap(yi)
    sigma_model_new = sigma_model_current(yi) + dsigma_y       

    call calc_pdf(sigma_model_current(yi), sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p0_sigma_y)
    call calc_pdf(sigma_model_new, sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p1_sigma_y)
	
    do jj=1,dim2   
        y_error_new(R_indices(:,jj)) = sigma_model_current(jj)    ! Provided dim2 isn't too big then should be fine
    enddo  

    ! Change one element of sigma_y
    y_error_new(R_indices(:,yi)) = sigma_model_new  ! R_indices = array of indices to which each sigma_y applies 
   			   
    sigma_y_new=sqrt(y_error_new**2 + sigma_measure**2)   

    n1T=sum((n0/sigma_y_new)**2)

    detval_new = sum(alog(sigma_y_new))           ! This is actually the inverse of the determinant

    ! Compute P1/P0 	
    pT= p1_sigma_y-p0_sigma_y - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant
    if (sigma_model_pdf .eq. 1) then
       if (sigma_model_new .lt. sigma_model_hparam1(yi)) pT = -1.e20
       if (sigma_model_new .gt. sigma_model_hparam2(yi)) pT = -1.e20
    endif
  
    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_model_current(yi) = sigma_model_new
       sigma_y_current = sigma_y_new
       detval_current = detval_new
       n0T=n1T
       !if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
       if(it .gt. burn_in) accept=accept + 1
       if(beta .eq. 1. .and. it .le. burn_in) accept_batch(yi)=accept_batch(yi) + 1
    else
       !;REJECT					
       !if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
       if(it .gt. burn_in) reject=reject + 1
       if(beta .eq. 1. .and. it .le. burn_in) reject_batch(yi)=reject_batch(yi) + 1
    endif

    if(beta .eq. 1 .and. it .le. burn_in .and. modulo(it,500) .eq. 1) then
       if (it .gt. 1) then
    !if(beta .eq. 1 .and. it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100) then
        av_acc =  max(0.01,1.0/sqrt(real(it)/500.))
        do ii=1,dim2 
          if(accept_batch(ii)+reject_batch(ii) .gt. 0) then          
             accep_prob = real(accept_batch(ii))/(accept_batch(ii)+reject_batch(ii))
             !if (ii .eq. 1) write(*,*) accep_prob
             if(accep_prob .lt. 0.2) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) - av_acc)
             if(accep_prob .gt. 0.6) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) + av_acc)
             accept_batch(ii) = 0
             reject_batch(ii) = 0
          endif
        enddo
       endif   ! it .gt. 1
    endif

n0T_out=n0T
sigma_y_out=sigma_y_current
sigma_model_out=sigma_model_current
accept_out=accept
reject_out=reject
detval_out=detval_current
stepsize_sig_out=stepsize_sigma_y
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE sigma_y_update




FUNCTION random_normal() RESULT(fn_val)
IMPLICIT NONE
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, half =0.5, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

SUBROUTINE closest_grid(region_v,lons,lats, plon,plat, np, ngrid)
implicit none

integer :: region_v(ngrid)
real    :: lons(ngrid), lats(ngrid), plon(np), plat(np)
integer :: lai, idx,pi
integer :: ngrid,np
real    :: maxdist, dist
    lai=1
    do idx =1,ngrid
   
        maxdist=1.e6
        do pi=1,np 
            dist=(lats(idx) - plat(pi))*(lats(idx) - plat(pi)) + (lons(idx) - plon(pi))*(lons(idx) - plon(pi))  
            ! Can't work out which way round??????????????
            if (dist .LT. maxdist) THEN
                region_v(lai)=pi
                maxdist=dist
            endif
        enddo
           
        lai=lai+1
    enddo

END SUBROUTINE closest_grid


SUBROUTINE closest_grid_sng(region,lon,lat, plon,plat, np)
implicit none

  integer :: region, pi, np
  real    :: plon(np), plat(np)
  real    :: lon, lat, maxdist, dist

  maxdist=1.e6
  do pi=1,np 
     dist=(lat - plat(pi))*(lat - plat(pi)) + (lon - plon(pi))*(lon - plon(pi))
     if (dist .LT. maxdist) THEN  
        region=pi
        maxdist=dist
     endif
  enddo

END SUBROUTINE closest_grid_sng

subroutine calc_pdf(x,pdf_param1,pdf_param2,pdf,p1)     ! Subroutine for calculating P1 for specified PDF

 Implicit none
 ! Args
 integer,intent(in)   :: pdf
 real,intent(in)      :: x, pdf_param1, pdf_param2
 real,intent(out)     :: p1
 integer              :: positive
 real                 :: mu, sigma
 ! Locals
real,parameter        :: pi = 3.14159265        ! a fixed constant

!! ALL PDFS ARE IN LOG-SPACE
!  PDF .eq. 'UNIFORM'
  If(pdf .eq. 1) then                 ! Uniform PDF = 1	

     if (x .lt. pdf_param1 .OR. x .gt. pdf_param2) then
        p1 = -1.e20
     else
        p1= alog(1./(pdf_param2 - pdf_param1))
     endif
     positive=0
  Endif

!  PDF .eq. 'GAUSSIAN'
  If(pdf .eq. 2) then                ! Gaussian PDF = 2
     p1=alog(1./(pdf_param2*sqrt(2.*pi)))+(-1.*(x - pdf_param1)**2 / 2./(pdf_param2**2))
     positive=0
  Endif

!  PDF .eq. 'LOGNORMAL'
  If(pdf .eq. 3) then                    ! Lognormal PDF = 3

    mu = alog(pdf_param1) - 0.5*alog(1. + pdf_param2**2/pdf_param1**2)
    sigma=sqrt(alog((pdf_param2/pdf_param1)**2 + 1.))


     p1=alog(1./(x*sigma*sqrt(2.*pi)))+( -1.*(alog(x) - mu)**2 / 2./(sigma**2))
     positive=1
  Endif

!  PDF .eq. 'EXPONENTIAL'
  If(pdf .eq. 4) then                    ! Exponential PDF = 4
     p1=alog(pdf_param1)+(-1.*pdf_param1*x)
     positive=1
  Endif


  if((x .le. 0.) .and. (positive .eq. 1)) p1 = -1.e20       ! Since log(0) isn't defined need to make p1 a very large negative number to ensure proposal is rejected

end subroutine calc_pdf

subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s, getpid
            integer(8) :: count, tms
           
            call random_seed(size = n)
            allocate(seed(n))

            ! First try if the OS provides a random number generator
            open(unit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
end subroutine init_random_seed

End Module mcmc_uncorr

