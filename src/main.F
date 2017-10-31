
!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE INITIAL_GRID(ng)
     USE GLOBAL  
     USE LOCAL   
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng
     INTEGER :: mbox1,nbox1,isk,mbox0,nbox0,m_move,n_move

     REAL(SP) :: Dphi_r,Dtheta_r  !ykchoi

# if defined (PARALLEL)
     INTEGER :: parent_area_i, parent_area_j
     INTEGER :: area_i, area_j
     INTEGER :: II,JJ
# endif
     
! for ng=1
     IF(ng==1)THEN
# if defined(CARTESIAN)
	 DX_Tile(ng) = DX
	 DY_Tile(ng) = DY
# else
       DX_Tile(1:Mloc,1:Nloc,ng) = DX(1:Mloc,1:Nloc)
       DY_Tile(1:Mloc,1:Nloc,ng) = DY(1:Mloc,1:Nloc)
# endif
       Depth_Tile(1:Mloc,1:Nloc,ng) = Depth(1:Mloc,1:Nloc)
       DepthX_Tile(1:Mloc1,1:Nloc,ng) = DepthX(1:Mloc1,1:Nloc)
       DepthY_Tile(1:Mloc,1:Nloc1,ng) = DepthY(1:Mloc,1:Nloc1)

       MASK_Tile(1:Mloc,1:Nloc,ng) = MASK(1:Mloc,1:Nloc)
       MASK9_Tile(1:Mloc,1:Nloc,ng) = MASK9(1:Mloc,1:Nloc)
       U_Tile(1:Mloc,1:Nloc,ng) = U(1:Mloc,1:Nloc)
       V_Tile(1:Mloc,1:Nloc,ng) = V(1:Mloc,1:Nloc)
       Ubar_Tile(1:Mloc,1:Nloc,ng) = Ubar(1:Mloc,1:Nloc)
       Vbar_Tile(1:Mloc,1:Nloc,ng) = Vbar(1:Mloc,1:Nloc)
       Eta_Tile(1:Mloc,1:Nloc,ng) = Eta(1:Mloc,1:Nloc)
# if defined(CARTESIAN)
!      We don't need latitude
# else
       Lat_Theta_Tile(1:Mloc,1:Nloc,ng) = Lat_Theta(1:Mloc,1:Nloc)
# endif

     ELSE

! for others
! mbox is box in a coarser grid, include ghostcells, m represents m+nghost

       !========================
	 ! case : ng >= 2
	 !========================

       isk=RATIO_SPACING(ng)
	 
       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

# if defined (PARALLEL)
	 EastParentID(:,:,ng)=-99;   WestParentID(:,:,ng)=-99;
	 SouthParentID(:,:,ng)=-99;  NorthParentID(:,:,ng)=-99;

	 parent_area_i = GridDimX( ng-1 )/px;
	 parent_area_j = GridDimY( ng-1 )/py;

	 if( n_west .eq. MPI_PROC_NULL ) THEN
	   DO J=1,Nloc
	      DO I=1,Nghost
               II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk !actually npx=0 at west
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk
			 !--- ykchoi(22/June/2017)
			 !These routines (II,JJ) are also used in subourtine INTERP_BC (bc_nesting.F).
		     !Thus, in ProcessorID, i,j are computed by using these II, JJ.
	         !Please see subourtine INTERP_BC before modifying these routines.

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )
	         
			 WestParentID(I,J,ng) = ProcessorID(area_i,area_j)
	      ENDDO
	   ENDDO
	 endif

	 if( n_east .eq. MPI_PROC_NULL ) THEN
	   DO J=1,Nloc
   	      DO I=Mloc-Nghost+1,Mloc
		     II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )

			 EastParentID(I+Nghost-Mloc,J,ng) = ProcessorID(area_i,area_j)  !be careful
	      ENDDO
	   ENDDO
	 endif

	 if( n_suth .eq. MPI_PROC_NULL ) then
	   DO I=1,Mloc
	      DO J=1,Nghost
               II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk !actually npy=0 at south

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )
	         
			 SouthParentID(I,J,ng) = ProcessorID(area_i,area_j)
	      ENDDO
	   ENDDO
	 endif

	 if( n_nrth .eq. MPI_PROC_NULL ) then
	   DO I=1,Mloc
	      DO J=Nloc-Nghost+1,Nloc
		     II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )

	         NorthParentID(I,J+Nghost-Nloc,ng) = ProcessorID(area_i,area_j)  !be careful
	      ENDDO
	   ENDDO
	 endif

# endif

! static
       !Careful when "Mloc" on grid1 < "Mloc" on grid2
	 ! --> Dx_Tile in some index = 0
       !DX_Tile(1:Mloc,1:Nloc,ng)=DX_Tile(1:Mloc,1:Nloc,ng-1)/isk
       !DY_Tile(1:Mloc,1:Nloc,ng)=DY_Tile(1:Mloc,1:Nloc,ng-1)/isk

# if defined (CARTESIAN)
       DX_Tile(ng) = DX_Tile(ng-1)/RATIO_SPACING(ng)
	 DY_Tile(ng) = DY_Tile(ng-1)/RATIO_SPACING(ng)
# else
       Dphi_r=Dphi*pi/180.0_SP
       Dtheta_r=Dtheta*pi/180.0_SP

       Do J=1,Nloc
       Do I=1,Mloc
# if defined (PARALLEL)
          Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                         +(npy*Nglob/py+J-1)*Dtheta_r
# else
          Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                        +(J-1)*Dtheta_r
# endif
          !DX_Tile(I,J,ng) = R_earth*Dphi_r*COS(Lat_theta(I,J))/isk
          !DY_Tile(I,J,ng) = R_earth*Dtheta_r/isk
	    DX_Tile(I,J,ng) = R_earth*Dphi_r*COS(Lat_theta(I,J))/TOTALRATIO_SPACING(ng)
	    DY_Tile(I,J,ng) = R_earth*Dtheta_r/TOTALRATIO_SPACING(ng)
       ENDDO
       ENDDO
# endif

! depth
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  Depth_Tile(:,:,ng-1),Depth_Tile(1:Mloc,1:Nloc,ng),ng)

       CALL LINEAR_INTERP(MaxDimX1,MaxDimY,Mloc1,Nloc,mbox1,nbox1,isk, &
                  DepthX_Tile(:,:,ng-1),DepthX_Tile(1:Mloc1,1:Nloc,ng),ng)

       CALL LINEAR_INTERP(MaxDimX,MaxDimY1,Mloc,Nloc1,mbox1,nbox1,isk, &
                  DepthY_Tile(:,:,ng-1),DepthY_Tile(1:Mloc,1:Nloc1,ng),ng)

! Eta
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  Eta_Tile(:,:,ng-1),Eta_Tile(1:Mloc,1:Nloc,ng),ng)

! calculate MASK and MASK9 separately 
       DO J=1,Nloc
       DO I=1,Mloc
        IF( Eta_Tile(I,J,ng) < -DEPTH_Tile(I,J,ng) ) THEN
         MASK_Tile(I,J,ng)=0
         Eta_Tile(I,J,ng)=MinDepth-Depth_Tile(I,J,ng)
        ELSE
         MASK_Tile(I,J,ng)=1
        ENDIF
       ENDDO
       ENDDO

     DO J=Jbeg-1,Jend+1
     DO I=Ibeg-1,Iend+1
      MASK9_Tile(I,J,ng)=MASK_Tile(I,J,ng)*MASK_Tile(I-1,J,ng)*MASK_Tile(I+1,J,ng)  &
                *MASK_Tile(I+1,J+1,ng)*MASK_Tile(I,J+1,ng)*MASK_Tile(I-1,J+1,ng) &
                *MASK_Tile(I+1,J-1,ng)*MASK_Tile(I,J-1,ng)*MASK_Tile(I-1,J-1,ng) 
      IF(ABS(Eta_Tile(I,J,ng))/MAX(DEPTH_Tile(I,J,ng),MinDepthFrc)>SWE_ETA_DEP)THEN
       MASK9_Tile(I,J,ng)=ZERO
      ENDIF

     ENDDO
     ENDDO

! U
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  U_Tile(:,:,ng-1),U_Tile(1:Mloc,1:Nloc,ng),ng)

! V
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  V_Tile(:,:,ng-1),V_Tile(1:Mloc,1:Nloc,ng),ng)

! Ubar
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  Ubar_Tile(:,:,ng-1),Ubar_Tile(1:Mloc,1:Nloc,ng),ng)

! Vbar
       CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  Vbar_Tile(:,:,ng-1),Vbar_Tile(1:Mloc,1:Nloc,ng),ng)

! Lat_theta
# if defined (CARTESIAN)
       
# else
	 CALL LINEAR_INTERP(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1,isk, &
                  Lat_Theta_Tile(:,:,ng-1),Lat_Theta_Tile(1:Mloc,1:Nloc,ng),ng)
# endif
! exchange for parallel !!!

     ENDIF  ! end ng

END SUBROUTINE INITIAL_GRID

!===========================================================================
!===========================================================================

# if defined (PARALLEL)
SUBROUTINE FIND_PROCESSOR( parent_area_i, PX, II, &
	                      parent_area_j, PY, JJ, area_i, area_j )
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: parent_area_i, PX, II, parent_area_j, PY, JJ
      INTEGER,INTENT(OUT) :: area_i, area_j     
      INTEGER :: K

	DO K=1,PX
	   IF( parent_area_i*(K-1) < II .and. II <= parent_area_i*K ) THEN
	       area_i = K
	   ENDIF
	ENDDO
	DO K=1,PY
	   IF( parent_area_j*(K-1) < JJ .and. JJ <= parent_area_j*K ) THEN
	       area_j = K
	   ENDIF
	ENDDO	         
ENDSUBROUTINE FIND_PROCESSOR
# endif

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE LINEAR_INTERP(MaxM,MaxN,M,N,mb,nb,isk,Fin,Fout,ng)
      USE PARAM
!ykchoi
# if defined (PARALLEL)
	USE GLOBAL, ONLY : GridDimX,GridDimY,Nghost,px,py,npx,npy
# endif

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,mb,nb,isk,MaxM,MaxN,ng   !ykchoi
      INTEGER :: II,JJ
      REAL(SP) :: rII,rJJ      
      REAL(SP),DIMENSION(MaxM,MaxN),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Fout

!ykchoi
# if defined (PARALLEL)
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGrid1  ! global including ghost
      INTEGER :: mm1, nn1, Mloc_grid1, Nloc_grid1
# endif

!ykchoi
# if defined (PARALLEL)
      mm1 = GridDimX(ng-1) + 2*Nghost
      nn1 = GridDimY(ng-1) + 2*Nghost
      ALLOCATE( VarGrid1(mm1,nn1) )

      Mloc_grid1 = GridDimX(ng-1)/px + 2*Nghost
      Nloc_grid1 = GridDimY(ng-1)/py + 2*Nghost 
      CALL GATHER_GRID( VarGrid1, Fin(1:Mloc_grid1,1:Nloc_grid1), &
               Mloc_grid1, Nloc_grid1, mm1, nn1, Nghost)
# endif


# if defined (PARALLEL)
       DO J=1,N
        DO I=1,M
         II = mb + (I-1+(M-2*Nghost)*npx)/isk
         JJ = nb + (J-1+(N-2*Nghost)*npy)/isk
         rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*npx)/REAL(isk) - REAL(II)
         rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*npy)/REAL(isk) - REAL(JJ)
         Fout(I,J)  &
	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ ) + &
             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
        ENDDO
       ENDDO
	 DEALLOCATE( VarGrid1 )
# else
! not parallel

       DO J=1,N
       DO I=1,M
        II=mb+(I-1)/isk
        JJ=nb+(J-1)/isk
        rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
        rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
        Fout(I,J)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
       ENDDO
       ENDDO
# endif

END SUBROUTINE LINEAR_INTERP

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi

SUBROUTINE LOAD_DATA(ng)

     USE GLOBAL  
     USE LOCAL

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng

# if defined(CARTESIAN)
      ! No variables "Lat_Theta"
# else
     Lat_Theta(1:Mloc,1:Nloc) = Lat_Theta_Tile(1:Mloc,1:Nloc,ng)
# endif

     Eta(1:Mloc,1:Nloc) = Eta_Tile(1:Mloc,1:Nloc,ng)
     Ubar(1:Mloc,1:Nloc) = Ubar_Tile(1:Mloc,1:Nloc,ng)
     Vbar(1:Mloc,1:Nloc) = Vbar_Tile(1:Mloc,1:Nloc,ng)  

     U(1:Mloc,1:Nloc) = U_Tile(1:Mloc,1:Nloc,ng)
     V(1:Mloc,1:Nloc) = V_Tile(1:Mloc,1:Nloc,ng)
     MASK(1:Mloc,1:Nloc) = MASK_Tile(1:Mloc,1:Nloc,ng)
     MASK9(1:Mloc,1:Nloc) = MASK9_Tile(1:Mloc,1:Nloc,ng)
     Depth(1:Mloc,1:Nloc) = Depth_Tile(1:Mloc,1:Nloc,ng)
     DepthX(1:Mloc1,1:Nloc) = DepthX_Tile(1:Mloc1,1:Nloc,ng)
     DepthY(1:Mloc,1:Nloc1) = DepthY_Tile(1:Mloc,1:Nloc1,ng)

# if defined(CARTESIAN)
     DX = DX_Tile(ng)
     DY = DY_Tile(ng)
# else
     DX(1:Mloc,1:Nloc) = DX_Tile(1:Mloc,1:Nloc,ng)
     DY(1:Mloc,1:Nloc) = DY_Tile(1:Mloc,1:Nloc,ng)
# endif

END SUBROUTINE LOAD_DATA

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE CALC_GRID(ng,NestStep,NestTotal)

     USE GLOBAL  
     USE LOCAL   
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng,NestStep,NestTotal

     INTEGER::ISTAGE

! Consider later (ykchoi)
! # if defined (SAMPLES)
!     IF(WaveMaker(1:7)=='LEF_SOL')THEN
!        CALL SOLITARY_WAVE_LEFT_BOUNDARY
!     ENDIF  
! # endif 

   ! update three variables
     Eta0=Eta
     Ubar0=Ubar
     Vbar0=Vbar

     CALL UPDATE_MASK

     CALL EXCHANGE(ng)    !AMR routine
     IF(ng.GT.1)THEN
       CALL USE_NESTING_BC(NestStep,NestTotal,ng)
     ENDIF

# if defined (COUPLING)
     CALL OneWayCoupling
# endif

! calculate other variables for fluxes
     H=Eta*Gamma3+Depth
     HU=H*U*MASK
     HV=H*V*MASK

     IF(ng.EQ.1)THEN
# if defined(CARTESIAN)
       CALL ESTIMATE_DT(Mloc,Nloc,DX,DY,  &
	      U(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),H(1:Mloc,1:Nloc),MinDepthFrc,DT,CFL,TIME)
# else
       CALL ESTIMATE_DT(Mloc,Nloc,DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	      U(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),H(1:Mloc,1:Nloc),MinDepthFrc,DT,CFL,TIME)
# endif
     ENDIF

     U0=U
     V0=V

! # if defined (VESSEL)
!       CALL VESSEL_FORCING    !consider later
! # endif

! 3-ORDER RUNGE-KUTTA TIME STEPPING
     DO ISTAGE=1,3

        IF(DISPERSION)THEN
          CALL Cal_Dispersion(ng)    !AMR routine
        ENDIF
	  
	  CALL FLUXES(ng)
	  
! # if defined (WIND)
!        IF(WindForce)THEN
!          CALL WIND_EFFECT    !consider later
!        ENDIF
! # endif

	  CALL SourceTerms(ng)   ! put sourceterms after fluxes in order to get eta_t 

	  CALL ESTIMATE_HUV(ISTAGE, ng)    !AMR routine

        CALL WAVE_BREAKING

        CALL EXCHANGE(ng)    !AMR routine

# if defined (COUPLING)
        CALL OneWayCoupling
# endif

        IF(ng.GT.1)THEN
          CALL USE_NESTING_BC(NestStep,NestTotal,ng)
        ENDIF

        IF(ng.Eq.1)THEN
          IF(WaveMaker(1:3)=='ABS') THEN
            CALL ABSORBING_GENERATING_BC
          ENDIF

          IF(DIRECT_SPONGE)THEN
            CALL SPONGE_DAMPING
          ENDIF
        ENDIF

     ENDDO

!==========================
! save for next time step
!==========================
     MASK_TILE(1:Mloc,1:Nloc,ng) = MASK(1:Mloc,1:Nloc)
     MASK9_TILE(1:Mloc,1:Nloc,ng) = MASK9(1:Mloc,1:Nloc)

     U_TILE(1:Mloc,1:Nloc,ng) = U(1:Mloc,1:Nloc)
     V_TILE(1:Mloc,1:Nloc,ng) = V(1:Mloc,1:Nloc)
     Ubar_TILE(1:Mloc,1:Nloc,ng) = Ubar(1:Mloc,1:Nloc)
     Vbar_TILE(1:Mloc,1:Nloc,ng) = Vbar(1:Mloc,1:Nloc)
       
     Eta_TILE(1:Mloc,1:Nloc,ng) = Eta(1:Mloc,1:Nloc)

!     CALL MIXING_STUFF    !consider later

!  find maximum eta velocity
	IF(ng.eq.1)THEN 
       
	 IF (OUT_Hmax.OR.OUT_Hmin.OR.OUT_Umax.OR.OUT_MFmax.OR.OUT_VORmax)THEN
         CALL MAX_MIN_PROPERTY
       ENDIF

	 CALL OUTPUT   !STATISTICS, STATIONS
       
	 CALL CHECK_BLOWUP

	ENDIF

END SUBROUTINE CALC_GRID

!===========================================================================
!===========================================================================

SUBROUTINE INITIALIZE_VARIABLES

     USE GLOBAL
     IMPLICIT NONE

     U_NESTING_EAST_PAR = U_NESTING_EAST
     V_NESTING_EAST_PAR = V_NESTING_EAST
     Z_NESTING_EAST_PAR = Z_NESTING_EAST

     U_NESTING_WEST_PAR = U_NESTING_WEST
     V_NESTING_WEST_PAR = V_NESTING_WEST
     Z_NESTING_WEST_PAR = Z_NESTING_WEST

     U_NESTING_SOUTH_PAR = U_NESTING_SOUTH
     V_NESTING_SOUTH_PAR = V_NESTING_SOUTH
     Z_NESTING_SOUTH_PAR = Z_NESTING_SOUTH

     U_NESTING_NORTH_PAR = U_NESTING_NORTH
     V_NESTING_NORTH_PAR = V_NESTING_NORTH
     Z_NESTING_NORTH_PAR = Z_NESTING_NORTH

      !Initialize
     Eta = ZERO
     Ubar = ZERO
     Vbar = ZERO

     U = ZERO
     V = ZERO
     MASK = ZERO
     MASK9 = ZERO
	     
     Depth = ZERO
     DepthX = ZERO
     DepthY = ZERO
     DX = ZERO
     DY = ZERO

     U_NESTING_EAST = ZERO
     V_NESTING_EAST = ZERO
     Z_NESTING_EAST = ZERO

     U_NESTING_WEST = ZERO
     V_NESTING_WEST = ZERO
     Z_NESTING_WEST = ZERO

     U_NESTING_SOUTH = ZERO
     V_NESTING_SOUTH = ZERO
     Z_NESTING_SOUTH = ZERO

     U_NESTING_NORTH = ZERO
     V_NESTING_NORTH = ZERO
     Z_NESTING_NORTH = ZERO

END SUBROUTINE INITIALIZE_VARIABLES

!===========================================================================
!===========================================================================