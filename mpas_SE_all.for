C pgf90 -mcmodel=medium -o mpas_SE.exe -s mpas_SE.for
      parameter(NX=800,NY=800,NZ=55,MR=400,MT=361,MRSE=200)  !'CHANGE ME'
      parameter(rearth=6371000.,pi=3.141592654)
      parameter(omega=7.292*1.0e-5,GG=9.80665)
      parameter(P0=100000.,cp=1004.5,epsi=0.609)
      double precision TX2,TY2,HTX,HTY,VTX,VTY,STX,STY,TTX,TTY,ATX,ATY
c      double precision ax,bx,cx
     
      INTEGER ITC,JTC,FN
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: 
     &  theta,pre,den,uu,vv,uu2,vv2,
     &  ww,qv,qc,qr,qi,qs,qg,rubl,rvbl,rucu,rvcu,rthbl,
     &  rthcu,rthlw,rthsw,rthdia,
     &  urt,vrt,wrt,prt,drt,theta_m,
     &  vr2,vt2,
     &  vr21,vt21,wrt1,prt1,drt1,
     &  HDIRT,THETART,QVRT,RUBLRT,RVBLRT,RTHBLRT,RBLR2,RBLT2,
     &  HDIRT1,THETART1,QVRT1,RUBLRT1,RVBLRT1,RTHBLRT1,RBLR21,RBLT21,
     &  WK2D,WK2D1,WK2DSE,WK2T,WK2D2,FV1,FV2,FV3,FV4,HT1,HT2,HT3,HT4,HT5
      REAL, DIMENSION(:,:), ALLOCATABLE :: 
     &  xlat,xlon,xland,terrain,
     &  TDX,TDY,DX,DY,uu2_temp,vv2_temp,
     &  bur,bvr,bwr,bpr,bdr,fcy,frt,
     &  bur1,bvr1,bwr1,bpr1,bdr1,fcy1,frt1,
     &  BHDIR,BTHETAR,BQVR,BRUBLR,BRVBLR,BRTHBLR,
     &  BHDIR1,BTHETAR1,BQVR1,BRUBLR1,BRVBLR1,BRTHBLR1,
     &  PII,VPT,CHI,AVTSE,ATI,CC,CHIC,BTI1,BTI2,XII,ETA,CTI,CON1,CON2,
     &  FV,HT,DSEST,SEST1,SEST2,SEDR,SEDZ,FVV,HTDR,HTDZ,
     &  SEUU,SEWW,DSEDZ,DSEDR,RBDR,FV11,FV22,FV33,FV44,
     &  HT11,HT22,HT33,HT44,HT55,
     &  DHTDR,DHTDZ,FVtal,DCHICDZ,USE,WSE,SS,GAM,SERSO1,SERSO2,SERSO3,
     &  SERSO4,SERSO5,SERSO,TWT,DPIIDR,PVRT,PVRT1,PVRT2,NRO,AGV1,AGV2,
     &  PIIb,VPTB,DPIIDRB,AVTSEB,CCB,WDOT,UDOT,UDOTB,WDOTB,
     &  CHIB,CHIUD,CHIUDB,CHIWD,CHIWDB,DUDDT,DWDDT,
     &  VDOT,VDDZ,DVDDZ,VTT,VTTb,terRT,terRT1,WK2D_ter1,WK2D_ter2
      REAL, DIMENSION(:), ALLOCATABLE :: BterR
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: JSUMALL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     &  DGAM,DA11,DA12,DB11,DB12,DB21,DB22,DC11,DC12,
     &  DSERSO1,DSERSO2,DSERSO3,DSERSO4,DSERSO,DSEST1,DSEST2
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: WK3D1,WK3D2
      REAL, DIMENSION(:,:), ALLOCATABLE :: Ushear0,Vshear0
      real Z00(NZ),DZ(NZ),DD(MR),BFR(MR),A1P(NZ)
      real TDX1(NX),TDY1(NY)
      real Ushear(NZ),Vshear(NZ)
      real DR,DTSE,SSmax,LDA
      double precision AIP12,AIM12,CKP12,CKM12,DRD,DZD,BTI1D1,BTI1D2,
     &  WSOD,DSS
      double precision DZDP(NZ)
      real*8 XCS,YCS
      character pmfile*50,pgfile*50,asfile*50,filename0*50
      character,dimension(5) :: filename*50
      character(len=50), ALLOCATABLE :: pfile(:) 
      real SPV
      integer ILS
      namelist /INP/DT,PGFILE,
     &              CLON,CLAT,CLONINI,CLATINI,KTC1,KTC2,
     &              NF,IGWB,IBaSE,ICYLI,IABC,INTE,IDZCONST,
     &              ISHEAR,IASY
      
      
      ALLOCATE(theta(NX,NY,NZ),
     &  pre(NX,NY,NZ),den(NX,NY,NZ),uu(NX,NY,NZ),vv(NX,NY,NZ),
     &  ww(NX,NY,NZ),qv(NX,NY,NZ),qc(NX,NY,NZ),qr(NX,NY,NZ),
     &  qi(NX,NY,NZ),qs(NX,NY,NZ),qg(NX,NY,NZ),rubl(NX,NY,NZ),
     &  rvbl(NX,NY,NZ),rucu(NX,NY,NZ),rvcu(NX,NY,NZ),rthbl(NX,NY,NZ),
     &  rthcu(NX,NY,NZ),rthlw(NX,NY,NZ),rthsw(NX,NY,NZ),
     &  rthdia(NX,NY,NZ),uu2(NX,NY,NZ),vv2(NX,NY,NZ))
      ALLOCATE(urt(MR,MT,NZ),vrt(MR,MT,NZ),
     &  wrt(MR,MT,NZ),drt(MR,MT,NZ),prt(MR,MT,NZ),
     &  wrt1(MR,MT,NZ),drt1(MR,MT,NZ),prt1(MR,MT,NZ),
     &  vr2(MR,MT,NZ),vt2(MR,MT,NZ),
     &  vr21(MR,MT,NZ),vt21(MR,MT,NZ),
     &  HDIRT(MR,MT,NZ),THETART(MR,MT,NZ),
     &  HDIRT1(MR,MT,NZ),THETART1(MR,MT,NZ),
     &  QVRT(MR,MT,NZ),RUBLRT(MR,MT,NZ),RVBLRT(MR,MT,NZ),
     &  QVRT1(MR,MT,NZ),RUBLRT1(MR,MT,NZ),RVBLRT1(MR,MT,NZ),
     &  RTHBLRT(MR,MT,NZ),RBLR2(MR,MT,NZ),RBLT2(MR,MT,NZ),
     &  RTHBLRT1(MR,MT,NZ),RBLR21(MR,MT,NZ),RBLT21(MR,MT,NZ),
     &  FV1(MR,MT,NZ),FV2(MR,MT,NZ),FV3(MR,MT,NZ),FV4(MR,MT,NZ),
     &  HT1(MR,MT,NZ),HT2(MR,MT,NZ),HT3(MR,MT,NZ),HT4(MR,MT,NZ),
     &  HT5(MR,MT,NZ))
      ALLOCATE(xlat(NX,NY),xlon(NX,NY),xland(NX,NY),terrain(NX,NY),
     &  TDX(NX,NY),TDY(NX,NY),DX(NX,NY),DY(NX,NY),
     &  uu2_temp(NX,NY),vv2_temp(NX,NY))
      ALLOCATE(bur(MR,NZ),bvr(MR,NZ),bwr(MR,NZ),bdr(MR,NZ),
     &  bur1(MR,NZ),bvr1(MR,NZ),bwr1(MR,NZ),bdr1(MR,NZ),
     &  bpr(MR,NZ),fcy(nx,ny),frt(MR,MT),BHDIR(MR,NZ),
     &  bpr1(MR,NZ),fcy1(nx,ny),frt1(MR,MT),BHDIR1(MR,NZ),
     &  BTHETAR(MR,NZ),BQVR(MR,NZ),BRUBLR(MR,NZ),BRVBLR(MR,NZ),
     &  BTHETAR1(MR,NZ),BQVR1(MR,NZ),BRUBLR1(MR,NZ),BRVBLR1(MR,NZ),
     &  BRTHBLR(MR,NZ),
     &  BRTHBLR1(MR,NZ),
     &  PII(MR,NZ),VPT(MR,NZ),CHI(MR,NZ),
     &  AVTSE(MR,NZ),ATI(MR,NZ),CC(MR,NZ),CHIC(MR,NZ),BTI1(MR,NZ),
     &  BTI2(MR,NZ),
     &  XII(MR,NZ),ETA(MR,NZ),CTI(MR,NZ),CON1(MR,NZ),CON2(MR,NZ),
     &  FV(MR,NZ),HT(MR,NZ),DSEST(MR,NZ),SEST1(MR,NZ),SEST2(MR,NZ),
     &  SEDR(MR,NZ),SEDZ(MR,NZ),FVV(MR,NZ),HTDZ(MR,NZ),HTDR(MR,NZ),
     &  SEUU(MR,NZ),SEWW(MR,NZ),DSEDR(MR,NZ),DSEDZ(MR,NZ),
     &  RBDR(MR,NZ),FV11(MR,NZ),FV22(MR,NZ),FV33(MR,NZ),FV44(MR,NZ),
     &  HT11(MR,NZ),
     &  HT22(MR,NZ),HT33(MR,NZ),HT44(MR,NZ),HT55(MR,NZ),DHTDR(MR,NZ),
     &  DHTDZ(MR,NZ),
     &  FVtal(MR,NZ),DCHICDZ(MR,NZ),USE(MR,NZ),WSE(MR,NZ),
     &  SS(MR,NZ),GAM(MR,NZ),SERSO1(MR,NZ),SERSO2(MR,NZ),
     &  SERSO3(MR,NZ),SERSO4(MR,NZ),SERSO5(MR,NZ),SERSO(MR,NZ),
     &  TWT(MR,NZ),DPIIDR(MR,NZ),
     &  PVRT(MR,NZ),PVRT1(MR,NZ),PVRT2(MR,NZ),NRO(MR,NZ),AGV1(MR,NZ),
     &  AGV2(MR,NZ),
     &  PIIB(MR,NZ),VPTB(MR,NZ),DPIIDRB(MR,NZ),AVTSEB(MR,NZ),CCB(MR,NZ),
     &  UDOTB(MR,NZ),WDOTB(MR,NZ),CHIB(MR,NZ),CHIUD(MR,NZ),
     &  CHIWD(MR,NZ),CHIUDB(MR,NZ),CHIWDB(MR,NZ),
     &  WDOT(MR,NZ),UDOT(MR,NZ),DUDDT(MR,NZ),DWDDT(MR,NZ),
     &  VDOT(MR,NZ),VDDZ(MR,NZ),DVDDZ(MR,NZ),
     &  VTT(MR,NZ),VTTb(MR,NZ),Ushear0(MR,NZ),Vshear0(MR,NZ))
      ALLOCATE(DGAM(MR,NZ),DSERSO1(MR,NZ),DSERSO2(MR,NZ),
     &  DSERSO3(MR,NZ),DSERSO4(MR,NZ),DSERSO(MR,NZ),
     &  DSEST1(MR,NZ),DSEST2(MR,NZ),DA11(MR,NZ),DA12(MR,NZ),DB11(MR,NZ),
     &  DB12(MR,NZ),DB21(MR,NZ),DB22(MR,NZ),DC11(MR,NZ),DC12(MR,NZ))
      ALLOCATE(WK2D(MR,NZ,50),WK2D1(MR,NZ,30),WK2D2(MR,NZ,30),
     &  WK3D1(MR,MT,NZ,15),WK3D2(MR,MT,NZ,15),
     &  WK2DSE(MRSE,NZ,10))
      ALLOCATE(terRT(MR,MT),terRT1(MR,MT),
     &  WK2D_ter1(MR,MT),WK2D_ter2(MR,MT))
      ALLOCATE(BterR(MR))
      ALLOCATE(JSUMALL(MR,NZ))

      filename0='All'
      filename(1)=trim(filename0)//'_WholeTC'
      filename(2)=trim(filename0)//'_UL'
      filename(3)=trim(filename0)//'_UR'
      filename(4)=trim(filename0)//'_DR'
      filename(5)=trim(filename0)//'_DL'
      pmfile='mpas_SE_z.inp'
c  read parameters from pmfile
      OPEN (UNIT=98,FILE=PMFILE,FORM='FORMATTED',STATUS='OLD')
      READ(98,NML=INP)
      CLOSE (UNIT=98)
C
       print *,'reading parameters from ',pmfile
      write(6,NML=INP)

      ALLOCATE(pfile(NF))
      open(unit=7,FILE="file.namelist")
      read(7,*) (pfile(i),i=1,NF) 
      print*,(pfile(i),i=1,NF)
      close (unit=7)

      ITC=((CLON-CLONINI)/0.025+1)  !shsar
      JTC=((CLAT-CLATINI)/0.025+1)  !shsar
      print*,' ITC=',ITC,' JTC=',JTC
 
      itime = 0
      
c      DO K=1,NZ
c      DO J=1,MT
c      DO I=1,MR
c      RTHLW(I,J,K)= 0.
c      RTHSW(I,J,K)= 0.
c      ENDDO
c      ENDDO
c      ENDDO

      DO K=1,NZ
      DO I=1,MR
      DO L= 1,50
      WK2D(I,K,L)= 0.
      ENDDO
      DO L= 1,30
      WK2D1(I,K,L)= 0.
      WK2D2(I,K,L)= 0.
      ENDDO
      DO J=1,MT
      DO L= 1,15
      WK3D1(I,J,K,L)= 0.
      WK3D2(I,J,K,L)= 0.
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      DO L= 1,10
      WK2DSE(I,K,L)= 0.
      ENDDO
      ENDDO
      ENDDO
 88   print *,'input filename (or exit) = ?'
      if (itime .eq. NF) then
c
      GO TO 99
c
      ELSE
C
      itime = itime + 1
      print *,PFILE(itime)

c      open(unit=33,file=pfile(itime),form='unformatted',
c     &  status='unknown',access='stream')
      open(unit=33,file=pfile(itime),form='unformatted',
     &  status='unknown')
      print*,'Reading variables start'
      read(33) (z00(k),k=1,NZ)
      read(33) ((xlat(i,j),i=1,NX),j=1,NY)
      read(33) ((xlon(i,j),i=1,NX),j=1,NY)
      read(33) ((xland(i,j),i=1,NX),j=1,NY)
      read(33) ((terrain(i,j),i=1,NX),j=1,NY)
      read(33) (((pre(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((theta(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((den(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((uu(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((vv(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((ww(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qv(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qc(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qr(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qi(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qs(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((qg(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rubl(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rvbl(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rucu(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rvcu(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rthbl(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rthcu(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rthlw(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rthsw(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      read(33) (((rthdia(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
c      read(33) (((theta_m(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
      CLOSE (UNIT=33)
      print*,'z00',z00
c      print*,'vv',vv(33,33,33)
c      print*,'vv',vv(413,413,3)
c      print*,'ww',ww(413,413,3)
c      print*,'theta',theta(413,413,3)
c      print*,'rthdia',rthdia(413,413,3)
      
c      DO I=1,NX
c      DO J=1,NY
c      uu(I,J,1)=0.
c      vv(I,J,1)=0.
c      ww(I,J,1)=0.
c      ENDDO
c      ENDDO
C  perform potential vorticity budgets on spheric coordinates
C  first compute grid mesh interval
      DO K=1,NZ-1
      DZ(K)= Z00(K+1)-Z00(K)
      ENDDO
      DZ(NZ)= DZ(NZ-1)
c  need to consider the DX and DY changing with LAT and LON...
      DO I=1,NX
      DO J=1,NY
      PHI= PI*XLAT(I,J)/180.
      TDX(I,J)= 2.*PI*rearth*cos(PHI)*XLON(I,J)/360.
      TDY(I,J)= rearth*PHI 
      ENDDO
      ENDDO
      DO I=1,NX-1
      DO J=1,NY
      DX(I,J)= TDX(I+1,J)-TDX(I,J)
      ENDDO
      ENDDO
      DO J=1,NY
      DX(NX,J)= DX(NX-1,J)
      ENDDO
      DO J=1,NY-1
      DO I=1,NX
      DY(I,J)= TDY(I,J+1)-TDY(I,J)
      ENDDO
      ENDDO
      DO I=1,NX
      DY(I,NY)= DY(I,NY-1)
      ENDDO
      DO I=1,NX
      TDX1(I) = TDX(I,JTC)
      ENDDO
      DO J=1,NY
      TDY1(J) = TDY(ITC,J)
      ENDDO

      print *,'get asymmetric parts with the center (ITC,JTC)'
      XCS = dble(TDX1(ITC))
      YCS = dble(TDY1(JTC))
      PRINT *,'The guessed XCS,YCS=',XCS,YCS

      DR = abs(TDX1(2)-TDX1(1))  
      DO I=1,MR
      DD(I) = 2*PI*DR*(I-1)/(MT-1)
      ENDDO
      DO I=1,NX
      DO J=1,NY
      PHI =  PI*XLAT(I,J)/180.
      FCY(I,J) = 2.*OMEGA*SIN(PHI)
      ENDDO
      ENDDO
      PRE= PRE*100.  ! converted to Pa
      Print*,'DR',DR   

c------------------------------------------------------------------
      IF (itime.eq.1 .and. IASY.eq.1) THEN
      DO K=KTC1,KTC2
      uu2_temp(:,:) = uu(:,:,K)
      vv2_temp(:,:) = vv(:,:,K)
      CALL ASYMV(uu2_temp,vv2_temp,NX,NY,TDX1,TDY1,XCS,YCS,0,0.,1,
     &       terRT,terrain)
      uu2(:,:,K) = uu2_temp
      vv2(:,:,K) = vv2_temp
      ENDDO

      open(unit=81,file=trim(pgfile)//'_asWind.dat',form='unformatted',
     &  access='direct',status='unknown',recl=NX*NY*NZ*4)
      WRITE(81,rec=1) uu2
      WRITE(81,rec=2) vv2
      CLOSE(UNIT=82)
      ENDIF

c------------------------------------------------------------------
      IF (ICYLI.EQ.1) THEN
      call CYLIN(terrain(:,:),NX,NY,TDX1,TDY1,terRT(:,:),MR,MT,
     &           BterR(:),DR,XCS,YCS,-999.,0,0,terRT(:,:),0.,0.,
     &           0,terrain(:,:),0,XLAND(:,:),XLON(:,:),CLON)

      DO K=KTC1,KTC2 !1,NZ
      print*,'K',K

      KM1= MAX0(K-1,1)
      ILS=0
      IF (Z00(K).LE.500.) THEN
        ILS=1
      ENDIF

      call CYLIN(UU(:,:,K),NX,NY,TDX1,TDY1,URT(:,:,K),MR,MT,BUR(:,K),DR,
     &           XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'URT',URT(3,3,K)
      call CYLIN(VV(:,:,K),NX,NY,TDX1,TDY1,VRT(:,:,K),MR,MT,BVR(:,K),DR,
     &           XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'VRT',VRT(3,3,K)
      call CYLIN(WW(:,:,K),NX,NY,TDX1,TDY1,WRT(:,:,K),MR,MT,BWR(:,K),DR,
     &           XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'WRT',WRT(3,3,K)
      call CYLIN(PRE(:,:,K),NX,NY,TDX1,TDY1,PRT(:,:,K),MR,MT,BPR(:,K),DR
     &           ,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'PRT',PRT(3,3,K)
      call CYLIN(DEN(:,:,K),NX,NY,TDX1,TDY1,DRT(:,:,K),MR,MT,BDR(:,K),DR
     &           ,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'DRT',DRT(3,3,K)
      call CYLIN(rthdia(:,:,K),NX,NY,TDX1,TDY1,HDIRT(:,:,K),MR,MT
     &         ,BHDIR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &         Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'HDIRT',HDIRT(3,3,K)
      call CYLIN(theta(:,:,K),NX,NY,TDX1,TDY1,THETART(:,:,K),MR,MT
     &       ,BTHETAR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &       Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'THETART',THETART(3,3,K)
      call CYLIN(qv(:,:,K),NX,NY,TDX1,TDY1,QVRT(:,:,K),MR,MT
     &          ,BQVR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &          Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'QVRT',QVRT(3,3,K)
      call CYLIN(rubl(:,:,K),NX,NY,TDX1,TDY1,RUBLRT(:,:,K),MR,MT
     &        ,BRUBLR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &        Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'RUBLRT',RUBLRT(3,3,K)
      call CYLIN(rvbl(:,:,K),NX,NY,TDX1,TDY1,RVBLRT(:,:,K),MR,MT
     &       ,BRVBLR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &       Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'RVBLRT',RVBLRT(3,3,K)
      call CYLIN(rthbl(:,:,K),NX,NY,TDX1,TDY1,RTHBLRT(:,:,K),MR,MT
     &       ,BRTHBLR(:,K),DR,XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),
     &       Z00(KM1),1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)
      print*,'RTHBLRT',RTHBLRT(3,3,K)
c      call CYLIN(RUCU(:,:,K),NX,NY,TDX1,TDY1,UCRT(:,:,K),MR,MT,BUCR(:,K)
c     &           ,DR,XCS,YCS,1.,0)
c      call CYLIN(RVCU(:,:,K),NX,NY,TDX1,TDY1,VCRT(:,:,K),MR,MT,BVCR(:,K)
c     &           ,DR,XCS,YCS,1.,0)
      ENDDO
      call CYLIN(FCY(:,:),NX,NY,TDX1,TDY1,FRT(:,:),MR,MT,BFR(:),DR,
     &           XCS,YCS,-999.,0,1,terRT(:,:),Z00(K),Z00(KM1),
     &           1,terrain(:,:),ILS,XLAND(:,:),XLON(:,:),CLON)

      call SVRVT(URT,VRT,VR2,VT2,MR,MT,NZ)
      call SVRVT(RUBLRT,RVBLRT,RBLR2,RBLT2,MR,MT,NZ)

      print*,'RTHBLRT',RTHBLRT(3,3,3)
      DO K=1,NZ
      DO I=2,MR
      TVR= 0.
      TVT= 0.
      JSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    VR2(I,J,K).NE.-999..and.VR2(IM1,J,K).NE.-999..and.
     &    VR2(IP1,J,K).NE.-999..and.VR2(I,JM1,K).NE.-999..and.
     &    VR2(I,JP1,K).NE.-999..and.VR2(I,J,KM1).NE.-999..and.
     &    VR2(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1
      TVR= TVR+ VR2(I,J,K)
      TVT= TVT+ VT2(I,J,K)
      ENDIF
      ENDDO
      BUR(I,K)= TVR/FLOAT(JSUM)
      BVR(I,K)= TVT/FLOAT(JSUM)
      ENDDO
      BUR(1,K)= VR2(1,1,K)
      BVR(1,K)= VT2(1,1,K)
      ENDDO
c
c chi,202312: calculate vertical wind shear
      DO K=1,NZ
      DO I=2,51         !shearEDIT
      TUshear= 0.
      TVshear= 0.
      IJSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    URT(I,J,K).NE.-999..and.URT(IM1,J,K).NE.-999..and.
     &    URT(IP1,J,K).NE.-999..and.URT(I,JM1,K).NE.-999..and.
     &    URT(I,JP1,K).NE.-999..and.URT(I,J,KM1).NE.-999..and.
     &    URT(I,J,KP1).NE.-999.) THEN
      IJSUM= IJSUM+1
      TUshear= TUshear+ URT(I,J,K)
      TVshear= TVshear+ VRT(I,J,K)
      ENDIF
      ENDDO
      Ushear0(I,K)= TUshear/FLOAT(IJSUM)
      Vshear0(I,K)= TVshear/FLOAT(IJSUM)
      ENDDO
      ENDDO

      DO K=1,NZ
      IJSUM= 0
      TUshear= 0.
      TVshear= 0.
      DO I=2,51
      IJSUM= IJSUM+(I-1)
      TUshear= TUshear+ Ushear0(I,K)*(I-1)
      TVshear= TVshear+ Vshear0(I,K)*(I-1)
      ENDDO
      Ushear(K)= TUshear/FLOAT(IJSUM)
      Vshear(K)= TVshear/FLOAT(IJSUM)
      ENDDO

c----------------------------------------------
      IF (itime .EQ. 1) THEN
      DO K=1,NZ
      DO I=1,MR
      WK2D1(I,K,1)= WK2D1(I,K,1) + BUR(I,K)
      WK2D1(I,K,2)= WK2D1(I,K,2) + BVR(I,K)
      WK2D1(I,K,3)= WK2D1(I,K,3) + BWR(I,K)
      WK2D1(I,K,4)= WK2D1(I,K,4) + BPR(I,K)
      WK2D1(I,K,5)= WK2D1(I,K,5) + BDR(I,K)
      WK2D1(I,K,6)= WK2D1(I,K,6) + BHDIR(I,K)
      WK2D1(I,K,7)= WK2D1(I,K,7) + BTHETAR(I,K)
      WK2D1(I,K,8)= WK2D1(I,K,8) + BQVR(I,K)
      WK2D1(I,K,9)= WK2D1(I,K,9) + FRT(I,K)
      WK2D1(I,K,10)= WK2D1(I,K,10) + BRUBLR(I,K)
      WK2D1(I,K,11)= WK2D1(I,K,11) + BRVBLR(I,K)
      WK2D1(I,K,12)= WK2D1(I,K,12) + BRTHBLR(I,K)
      DO J=1,MT
      WK3D1(I,J,K,1)= WK3D1(I,J,K,1) + VR2(I,J,K)   
      WK3D1(I,J,K,2)= WK3D1(I,J,K,2) + VT2(I,J,K)
      WK3D1(I,J,K,3)= WK3D1(I,J,K,3) + WRT(I,J,K)
      WK3D1(I,J,K,4)= WK3D1(I,J,K,4) + PRT(I,J,K)
      WK3D1(I,J,K,5)= WK3D1(I,J,K,5) + DRT(I,J,K)
      WK3D1(I,J,K,6)= WK3D1(I,J,K,6) + HDIRT(I,J,K)
      WK3D1(I,J,K,7)= WK3D1(I,J,K,7) + THETART(I,J,K)
      WK3D1(I,J,K,8)= WK3D1(I,J,K,8) + QVRT(I,J,K)
      WK3D1(I,J,K,9)= WK3D1(I,J,K,9) + RBLR2(I,J,K)
      WK3D1(I,J,K,10)= WK3D1(I,J,K,10) + RBLT2(I,J,K)
      WK3D1(I,J,K,11)= WK3D1(I,J,K,11) + RTHBLRT(I,J,K)
      WK2D_ter1(I,J) = terRT(I,J) 
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO K=1,NZ
      DO I=1,MR
      WK2D2(I,K,1)= WK2D2(I,K,1) + BUR(I,K)
      WK2D2(I,K,2)= WK2D2(I,K,2) + BVR(I,K)
      WK2D2(I,K,3)= WK2D2(I,K,3) + BWR(I,K)
      WK2D2(I,K,4)= WK2D2(I,K,4) + BPR(I,K)
      WK2D2(I,K,5)= WK2D2(I,K,5) + BDR(I,K)
      WK2D2(I,K,6)= WK2D2(I,K,6) + BHDIR(I,K)
      WK2D2(I,K,7)= WK2D2(I,K,7) + BTHETAR(I,K)
      WK2D2(I,K,8)= WK2D2(I,K,8) + BQVR(I,K)
      WK2D2(I,K,9)= WK2D2(I,K,9) + FRT(I,K)
      WK2D2(I,K,10)= WK2D2(I,K,10) + BRUBLR(I,K)
      WK2D2(I,K,11)= WK2D2(I,K,11) + BRVBLR(I,K)
      WK2D2(I,K,12)= WK2D2(I,K,12) + BRTHBLR(I,K)
      DO J=1,MT
      WK3D2(I,J,K,1)= WK3D2(I,J,K,1) + VR2(I,J,K)
      WK3D2(I,J,K,2)= WK3D2(I,J,K,2) + VT2(I,J,K)
      WK3D2(I,J,K,3)= WK3D2(I,J,K,3) + WRT(I,J,K)
      WK3D2(I,J,K,4)= WK3D2(I,J,K,4) + PRT(I,J,K)
      WK3D2(I,J,K,5)= WK3D2(I,J,K,5) + DRT(I,J,K)
      WK3D2(I,J,K,6)= WK3D2(I,J,K,6) + HDIRT(I,J,K)
      WK3D2(I,J,K,7)= WK3D2(I,J,K,7) + THETART(I,J,K)
      WK3D2(I,J,K,8)= WK3D2(I,J,K,8) + QVRT(I,J,K)
      WK3D2(I,J,K,9)= WK3D2(I,J,K,9) + RBLR2(I,J,K)
      WK3D2(I,J,K,10)= WK3D2(I,J,K,10) + RBLT2(I,J,K)
      WK3D2(I,J,K,11)= WK3D2(I,J,K,11) + RTHBLRT(I,J,K)
      WK2D_ter2(I,J) = terRT(I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      print*,'BUR',BUR(50,50)
     
      open(unit=11,file=trim(pgfile)//'_cyli1_1.dat',form='unformatted',
     &  access='direct',status='unknown',recl=MR*NZ*4)
      DO L=1,12
      WRITE(11,rec=L) WK2D1(:,:,L)
      ENDDO
      CLOSE (UNIT=11)
     
      open(unit=12,file=trim(pgfile)//'_cyli1_2.dat',form='unformatted',
     &  access='direct',status='unknown',recl=MR*MT*NZ*4)
      DO L=1,11
      WRITE(12,rec=L) WK3D1(:,:,:,L)
      ENDDO
      CLOSE (UNIT=12)
    
      open(unit=13,file=trim(pgfile)//'_cyli1_z.dat',form='unformatted',
     &  access='direct',status='unknown',recl=NZ*4)
      WRITE(13,rec=1) Z00
      CLOSE (UNIT=13)

      open(unit=14,file=trim(pgfile)//'_cyli1_ter.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MR*MT*4)
      WRITE(14,rec=1) WK2D_ter1(:,:)
      CLOSE (UNIT=14)

      open(unit=15,file=trim(pgfile)//'_cyli2_1.dat',form='unformatted',
     &  access='direct',status='unknown',recl=MR*NZ*4)
      DO L=1,12
      WRITE(15,rec=L) WK2D2(:,:,L)
      ENDDO
      CLOSE (UNIT=15)

      open(unit=16,file=trim(pgfile)//'_cyli2_2.dat',form='unformatted',
     &  access='direct',status='unknown',recl=MR*MT*NZ*4)
      DO L=1,11
      WRITE(16,rec=L) WK3D2(:,:,:,L)
      ENDDO
      CLOSE (UNIT=16)
      
      open(unit=17,file=trim(pgfile)//'_cyli2_z.dat',form='unformatted',
     &  access='direct',status='unknown',recl=NZ*4)
      WRITE(17,rec=1) Z00
      CLOSE (UNIT=17)

      open(unit=18,file=trim(pgfile)//'_cyli2_ter.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MR*MT*4)
      WRITE(18,rec=1) WK2D_ter2(:,:)
      CLOSE (UNIT=18)

      open(unit=19,file=trim(pgfile)//'_cyli2_shear.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=NZ*4)
      WRITE(19,rec=1) Ushear
      WRITE(19,rec=2) Vshear
      CLOSE (UNIT=19)

      ENDIF
c
      GO TO 88
c
      ENDIF
c
  99  CONTINUE

c------------------------------------------------------------------
      print*,'ABC',ABC
      IF (IABC.EQ.1) THEN
      open(unit=22,file=trim(pgfile)//'_cyli2_1.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(22) ((BUR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BVR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BWR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BPR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BDR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BHDIR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BTHETAR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BQVR(i,k),i=1,MR),k=1,NZ)
         read(22) ((FRT(i,k),i=1,MR),k=1,NZ)
         read(22) ((BRUBLR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BRVBLR(i,k),i=1,MR),k=1,NZ)
         read(22) ((BRTHBLR(i,k),i=1,MR),k=1,NZ)
      CLOSE (UNIT=22)
      print*,'BUR33',BUR(3,3)

      open(unit=23,file=trim(pgfile)//'_cyli2_2.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(23) (((VR2(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((VT2(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((WRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((PRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((DRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((HDIRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((THETART(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((QVRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((RBLR2(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((RBLT2(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(23) (((RTHBLRT(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
      CLOSE (UNIT=23)

      open(unit=24,file=trim(pgfile)//'_cyli1_1.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(24) ((BUR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BVR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BWR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BPR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BDR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BHDIR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BTHETAR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BQVR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((FRT1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BRUBLR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BRVBLR1(i,k),i=1,MR),k=1,NZ)
         read(24) ((BRTHBLR1(i,k),i=1,MR),k=1,NZ)
      CLOSE (UNIT=24)
      print*,'BUR1',BUR1(3,3)
      
      open(unit=25,file=trim(pgfile)//'_cyli1_2.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(25) (((VR21(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((VT21(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((WRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((PRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((DRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((HDIRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((THETART1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((QVRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((RBLR21(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((RBLT21(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
         read(25) (((RTHBLRT1(i,j,k),i=1,MR),j=1,MT),k=1,NZ)
      CLOSE (UNIT=25)
      print*,'VT2',VT21(3,3,3)

      open(unit=26,file=trim(pgfile)//'_cyli1_z.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(26) (Z00(k),k=1,NZ)
      CLOSE (UNIT=26)

      open(unit=27,file=trim(pgfile)//'_cyli2_ter.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(27) ((terRT(i,j),i=1,MR),j=1,MT)
      CLOSE (UNIT=27)

      open(unit=28,file=trim(pgfile)//'_cyli1_ter.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(28) ((terRT1(i,j),i=1,MR),j=1,MT)
      CLOSE (UNIT=28)

      open(unit=29,file=trim(pgfile)//'_cyli2_shear.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(29) (Ushear(k),k=1,NZ)
         read(29) (Vshear(k),k=1,NZ)
      CLOSE (UNIT=29)
      
c---------------------------------------------------------------
c define the Exner pressure and virtual potential temperature
      DO K=1,NZ
      DO I=1,MRSE
      PII2= 0.
      VPT2= 0.
      JSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    PRT(I,J,K).NE.-999..and.PRT(IM1,J,K).NE.-999..and.
     &    PRT(IP1,J,K).NE.-999..and.PRT(I,JM1,K).NE.-999..and.
     &    PRT(I,JP1,K).NE.-999..and.PRT(I,J,KM1).NE.-999..and.
     &    PRT(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1
c     Exner pressure
      PII1= cp*(PRT(I,J,K)/P0)**0.286
      PII2= PII2 + PII1
c     virtual potential temperature
      VPT1= THETART(I,J,K)*(1+epsi*QVRT(I,J,K))
      VPT2= VPT2 + VPT1
      ENDIF
      ENDDO
      PII(I,K)= PII2/FLOAT(JSUM)
      VPT(I,K)= VPT2/FLOAT(JSUM)
      CHI(I,K)= 1./VPT(I,K)
      ENDDO
      ENDDO
            
      PRINT*,'PII',PII(3,1)
      PRINT*,'VPT',VPT(3,1)
      
      DO K=1,NZ
      DO I=1,MRSE
      VT22= 0.
      JSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    VT2(I,J,K).NE.-999..and.VT2(IM1,J,K).NE.-999..and.
     &    VT2(IP1,J,K).NE.-999..and.VT2(I,JM1,K).NE.-999..and.
     &    VT2(I,JP1,K).NE.-999..and.VT2(I,J,KM1).NE.-999..and.
     &    VT2(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1
      VT23= VT2(I,J,K)
      VT22= VT22 + VT23
      ENDIF
      ENDDO
      VTT(I,K)= VT22/FLOAT(JSUM)
      ENDDO
      ENDDO
      DO I=1,MRSE
      VTT(I,1)= 2*VTT(I,2)-VTT(I,3)
      ENDDO
c--------------------------------------------------------
c gradient-wind balance
      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR = FLOAT(I-1)*DR
      AVT1 = -(FRT(I,1)*RR)/2
      AVT2 = AVT1*AVT1
      DPIIDR(I,K) = FNDIF1S(PII(IM1,K),PII(I,K),PII(IP1,K),DR1,DR2,SR)
      AVT3 = RR*VPT(I,K)*DPIIDR(I,K)
      IF (DPIIDR(I,K) .GT.0) THEN
      AVTSE(I,K) = AVT1+SQRT(AVT2+AVT3)
      ELSE
      AVTSE(I,K) = AVT1+SQRT(AMAX1(0.,AVT2+AVT3))
      ENDIF
      AVTSE(I,K) = AMAX1(0.,AVTSE(I,K))

c      IF (DPIIDR(I,K) .GT. 0) THEN
c      AVTSE(I,K)= AVT1+SQRT(AVT2+AVT3)
c      ELSE IF (DPIIDR(I,K) .EQ. 0) THEN
c      AVTSE(I,K)= 0.
c      ELSE 
c      AVTSE(I,K)= AVT1+SQRT(AMAX1(0.,AVT2+AVT3))
c      ENDIF
      ENDDO
      ENDDO

      DO K=1,NZ
      AVTSE(1,K)=0.
      ENDDO
      PRINT*,'AVTSE',AVTSE(3,1)
      
      IF (IGWB .NE. 1) THEN
      DO K=1,NZ
      DO I=1,MRSE
      AVTSE(I,K)= VTT(I,K)
      ENDDO
      ENDDO
      ENDIF
      
      DO K=1,NZ
      DO I=2,MRSE
      RR = FLOAT(I-1)*DR
      CC(I,K)=AVTSE(I,K)*AVTSE(I,K)/RR + AVTSE(I,K)*FRT(I,1)
      ENDDO
      ENDDO

      DO K=1,NZ
      CC(1,K)= 0.
      DPIIDR(1,K)= DPIIDR(2,K)
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      CHIC(I,K)=CC(I,K)*CHI(I,K)
      ENDDO
      ENDDO
c-----------------------------------------------------
c calculate A, B1, B2, and C
      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MR)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR = FLOAT(I-1)*DR
c A~
      DCHIDZ= FNDIF1S(CHI(I,KM1),CHI(I,K),CHI(I,KP1),DZ1,DZ2,SZ)
      NRO1= -GG/CHI(I,K)
      RN2= NRO1*DCHIDZ
      RN= AMAX1(1.0e-5,RN2)
      ATI(I,K)= RN*CHI(I,K)/(RR*BDR(I,K))
c B1~
      DCHIDR= FNDIF1S(CHI(IM1,K),CHI(I,K),CHI(IP1,K),DR1,DR2,SR)    
      BTI11= 1./(RR*BDR(I,K))
      BTI1(I,K)= GG*BTI11*DCHIDR
c B2~
      DCHICDZ(I,K)= FNDIF1S(CHIC(I,KM1),CHIC(I,K),CHIC(I,KP1),
     $  DZ1,DZ2,SZ)
      BTI2(I,K)= -BTI11*DCHICDZ(I,K)

c C~
      CTI1= 1./(RR*BDR(I,K))
      XII(I,K)= FRT(I,1)+2.*AVTSE(I,K)/RR
      DVDR= FNDIF1S(AVTSE(IM1,K),AVTSE(I,K),AVTSE(IP1,K),DR1,DR2,SR)
      ETA(I,K)= FRT(I,1)+DVDR+AVTSE(I,K)/RR
      CTI2= XII(I,K)*ETA(I,K)*CHI(I,K)
      DCHIDR= FNDIF1S(CHI(IM1,K),CHI(I,K),CHI(IP1,K),DR1,DR2,SR)
      CTI3= CC(I,K)*DCHIDR
      CTI(I,K)= CTI1*(CTI2+CTI3)
      IF (CTI(I,K) .LE. 0.) THEN
      CTI(I,K)= -0.001*CTI(I,K)
      ENDIF
c      CTI(I,K)= AMAX1(1.0e-12,abs(CTI(I,K)))
      ENDDO
      ENDDO
c
c deal with ATI near the inner boundary
      DO K=1,NZ
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DCHIDZ1= FNDIF1S(CHI(1,KM1),CHI(1,K),CHI(1,KP1),DZ1,DZ2,SZ)
      DCHIDZ2= FNDIF1S(CHI(2,KM1),CHI(2,K),CHI(2,KP1),DZ1,DZ2,SZ)
      DCHIDZ1P= (DCHIDZ1+DCHIDZ2)/2.
      BDR1P= (BDR(1,K)+BDR(2,K))/2.
      CHI1P= (CHI(1,K)+CHI(2,K))/2.
      RR1P = DR/2.
c      A1P(K)= -GG/(RR1P*BDR1P)*DCHIDZ1P
c      A1P(K)= 1.0e-4*CHI1P/(RR1P*BDR1P)
      
      NRO1= -GG/CHI1P
      RN2= NRO1*DCHIDZ1P
      RN= AMAX1(1.0e-5,RN2)
      A1P(K)= RN*CHI1P/(RR1P*BDR1P)
      ENDDO

       
      DO K=1,NZ
      ATI(1,K)= ATI(2,K)
      BTI1(1,K)= BTI1(2,K)
      BTI2(1,K)= BTI2(2,K)
      CTI(1,K)= CTI(2,K)
      ETA(1,K)= ETA(2,K)
      ENDDO

      IF (IBaSE .EQ. 1) THEN
      DO K=1,NZ
      DO I=1,MRSE
      BTI1(I,K)= BTI2(I,K)
      ENDDO
      ENDDO
      ENDIF
   
      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MR)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      NRO1= -GG/CHI(I,K)
      DCHIDZ= FNDIF1S(CHI(I,KM1),CHI(I,K),CHI(I,KP1),DZ1,DZ2,SZ)
      NRO(I,K)= SQRT(AMAX1(0.,NRO1*DCHIDZ))
      ENDDO
      ENDDO
      call MAXMIN2("NRO",NRO,MRSE,1,NZ,MR,NZ,0,0)

      PRINT*,'BTI1(3,3)',BTI1(3,3)
      call MAXMIN2("1:ATI",ATI,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("1:BTI1",BTI1,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("1:BTI2",BTI2,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("1:CTI",CTI,MRSE,1,NZ,MR,NZ,0,0)
c B^2-AC<0
      DO K=1,NZ
      DO I=1,MRSE
      CON1(I,K)= BTI1(I,K)*BTI1(I,K)-ATI(I,K)*CTI(I,K)
      IF (CON1(I,K) .GT. 0.) THEN
      BTI1(I,K)= 0.98*SQRT(ATI(I,K)*CTI(I,K))
c      BTI1(I,K)= 0.5*BTI1(I,K)
      ENDIF
      CON2(I,K)= BTI2(I,K)*BTI2(I,K)-ATI(I,K)*CTI(I,K)
      IF (CON2(I,K) .GT. 0.) THEN
      BTI2(I,K)= 0.98*SQRT(ATI(I,K)*CTI(I,K))
c      BTI2(I,K)= 0.5*BTI2(I,K)
      ENDIF
      CON1(I,K)= BTI1(I,K)*BTI1(I,K)-ATI(I,K)*CTI(I,K)
      CON2(I,K)= BTI2(I,K)*BTI2(I,K)-ATI(I,K)*CTI(I,K)
      ENDDO
      ENDDO

      call MAXMIN2("CON1",CON1,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("CON2",CON2,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("2:ATI",ATI,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("2:BTI1",BTI1,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("2:BTI2",BTI2,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("2:CTI",CTI,MRSE,1,NZ,MR,NZ,0,0)
c----------------------------------------------------
c calculate potential vorticity
      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MR)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      DVPTDZ= FNDIF1S(VPT(I,KM1),VPT(I,K),VPT(I,KP1),
     &  DZ1,DZ2,SZ)
      DVPTDR= FNDIF1S(VPT(IM1,K),VPT(I,K),VPT(IP1,K),
     &  DR1,DR2,SR)
      DVDR= FNDIF1S(AVTSE(IM1,K),AVTSE(I,K),AVTSE(IP1,K),DR1,DR2,SR)
      DVDZ= FNDIF1S(AVTSE(I,KM1),AVTSE(I,K),AVTSE(I,KP1),DZ1,DZ2,SZ)
      PVRT1(I,K)= ETA(I,K)*DVPTDZ
      PVRT2(I,K)= -DVDZ*DVPTDR
      PVRT(I,K)= (PVRT1(I,K)+PVRT2(I,K))/BDR(I,K)
      ENDDO
      ENDDO
      call MAXMIN2("PVRT",PVRT,MRSE,1,NZ,MR,NZ,1,0)      
    
c-------------------------------------------------------------------
c Udot and Wdot source
      DO K=1,NZ
      DO I=1,MRSE
      PII2= 0.
      VPT2= 0.
      JSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    PRT1(I,J,K).NE.-999..and.PRT1(IM1,J,K).NE.-999..and.
     &    PRT1(IP1,J,K).NE.-999..and.PRT1(I,JM1,K).NE.-999..and.
     &    PRT1(I,JP1,K).NE.-999..and.PRT1(I,J,KM1).NE.-999..and.
     &    PRT1(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1
c     Exner pressure
      PII1= cp*(PRT1(I,J,K)/P0)**0.286
      PII2= PII2 + PII1
c     virtual potential temperature
      VPT1= THETART1(I,J,K)*(1+epsi*QVRT1(I,J,K))
      VPT2= VPT2 + VPT1
      ENDIF
      ENDDO
      PIIb(I,K)= PII2/FLOAT(JSUM)
      VPTb(I,K)= VPT2/FLOAT(JSUM)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      VT22b= 0.
      JSUM= 0
      DO J=1,MT-1
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    VT21(I,J,K).NE.-999..and.VT21(IM1,J,K).NE.-999..and.
     &    VT21(IP1,J,K).NE.-999..and.VT21(I,JM1,K).NE.-999..and.
     &    VT21(I,JP1,K).NE.-999..and.VT21(I,J,KM1).NE.-999..and.
     &    VT21(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1
      VT23b= VT21(I,J,K)
      VT22b= VT22b + VT23b
      ENDIF
      ENDDO
      VTTb(I,K)= VT22b/FLOAT(JSUM)
      ENDDO
      ENDDO
      DO I=1,MRSE
      VTTb(I,1)= 2*VTTb(I,2)-VTTb(I,3)
      ENDDO      
      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR = FLOAT(I-1)*DR
c gradient-wind balance with the model output
      AVT1 = -(FRT1(I,1)*RR)/2
      AVT2 = AVT1*AVT1
      DPIIDRb(I,K) = FNDIF1S(PIIb(IM1,K),PIIb(I,K),PIIb(IP1,K),
     $  DR1,DR2,SR)
      AVT3 = RR*VPTb(I,K)*DPIIDRb(I,K)
      IF (DPIIDRb(I,K) .GT.0) THEN
      AVTSEb(I,K) = AVT1+SQRT(AVT2+AVT3)
      ELSE
      AVTSEb(I,K) = AVT1+SQRT(AMAX1(0.,AVT2+AVT3))
      ENDIF
      AVTSEb(I,K) = AMAX1(0.,AVTSEb(I,K))
C
c      IF (DPIIDRb(I,K) .GT. 0) THEN
c      AVTSEb(I,K)= AVT1+SQRT(AVT2+AVT3)
c      ELSE IF (DPIIDRb(I,K) .EQ. 0) THEN
c      AVTSEb(I,K)= 0.
c      ELSE
c      AVTSEb(I,K)= AVT1+SQRT(AMAX1(0.,AVT2+AVT3))
c      ENDIF
      ENDDO
      ENDDO

      DO K=1,NZ
      AVTSEb(1,K)=0.
      ENDDO

      IF (IGWB .NE. 1) THEN
      DO K=1,NZ
      DO I=1,MRSE
      AVTSEb(I,K)= VTTb(I,K)
      ENDDO
      ENDDO
      ENDIF

      PRINT*,'AVTSE',AVTSEb(3,1)
      DO K=1,NZ
      DO I=2,MRSE
      RR = FLOAT(I-1)*DR
      CCb(I,K)=AVTSEb(I,K)*AVTSEb(I,K)/RR + AVTSEb(I,K)*FRT1(I,1)
      ENDDO
      ENDDO

      DO K=1,NZ
      CCb(1,K)= 0.
      DPIIDRb(1,K)= DPIIDRb(2,K)
      ENDDO

c Udot
      DO K=1,NZ
      DO I=1,MRSE
      UDOT(I,K)= CC(I,K)-VPT(I,K)*DPIIDR(I,K)
      UDOTb(I,K)= CCb(I,K)-VPTb(I,K)*DPIIDRb(I,K)
      CHIb(I,K)= 1./VPTb(I,K)
      CHIUD(I,K)= CHI(I,K)*UDOT(I,K)
      CHIUDb(I,K)= CHIb(I,K)*UDOTb(I,K)
      ENDDO
      ENDDO

c Wdot
      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MR)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      DPIIDZ= FNDIF1S(PII(I,KM1),PII(I,K),PII(I,KP1),
     &  DZ1,DZ2,SZ)
      DPIIDZb= FNDIF1S(PIIb(I,KM1),PIIb(I,K),PIIb(I,KP1),
     &  DZ1,DZ2,SZ)
      WDOT(I,K)= -DPIIDZ*VPT(I,K)-GG
      WDOTb(I,K)= -DPIIDZb*VPTb(I,K)-GG
      CHIWD(I,K)= CHI(I,K)*WDOT(I,K)
      CHIWDb(I,K)= CHIb(I,K)*WDOTb(I,K)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MR)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      DCHIUDDZ= FNDIF1S(CHIUD(I,KM1),CHIUD(I,K),CHIUD(I,KP1),
     &  DZ1,DZ2,SZ)
      DCHIUDbDZ= FNDIF1S(CHIUDb(I,KM1),CHIUDb(I,K),CHIUDb(I,KP1),
     &  DZ1,DZ2,SZ)
      DCHIWDDR= FNDIF1S(CHIWD(IM1,K),CHIWD(I,K),CHIWD(IP1,K),
     &  DR1,DR2,SR)
      DCHIWDbDR= FNDIF1S(CHIWDb(IM1,K),CHIWDb(I,K),CHIWDb(IP1,K),
     &  DR1,DR2,SR)
      DUDDT(I,K)= (DCHIUDDZ-DCHIUDbDZ)/DT
      DWDDT(I,K)= (DCHIWDDR-DCHIWDbDR)/DT
CCC      DUDDT(I,K)= 0.
CCC      DWDDT(I,K)= 0.
      ENDDO
      ENDDO

c chi,202312: calculate shear and determine 4 quadrant
      FN= 1
      IF (ISHEAR.EQ.1) THEN
        UshearF= Ushear(29)- Ushear(17) !shearEDIT
        VshearF= Vshear(29)- Vshear(17)      
        RAtoDE= 180./PI 
        degree= RAtoDE*ACOS(UshearF/SQRT(UshearF**2.+VshearF**2.))
        IF (VshearF.LT.0.) degree= 360.- degree
        SJB= INT(degree)

        IQUAR= 0
  302   SJB= SJB+ 90
        IF (SJB.GT.360) SJB= SJB-360
        SJE= SJB+ 89
        FN= FN+ 1
        IQUAR= IQUAR+ 1
        IF (IQUAR.GT.4) GO TO 303
        print*,'quadrant:',SJB,SJE
        GO TO 301  
      ELSE
        SJB= 1
        SJE= MT-1
        GO TO 301
      ENDIF

c momentum source
  301 DO K=1,NZ
      DO I=2,MRSE
      JSUMALL(I,K)=0
      DO J=1,MT-1
      FV1(I,J,K)= 0.
      FV2(I,J,K)= 0.
      FV3(I,J,K)= 0.
      FV4(I,J,K)= 0.
      HT1(I,J,K)= 0.
      HT2(I,J,K)= 0.
      HT3(I,J,K)= 0.
      HT4(I,J,K)= 0.
      HT5(I,J,K)= 0.

      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)

      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    VT2(I,J,K).NE.-999..and.VT2(IM1,J,K).NE.-999..and.
     &    VT2(IP1,J,K).NE.-999..and.VT2(I,JM1,K).NE.-999..and.
     &    VT2(I,JP1,K).NE.-999..and.VT2(I,J,KM1).NE.-999..and.
     &    VT2(I,J,KP1).NE.-999.) THEN
      JSUMALL(I,K)= JSUMALL(I,K)+1
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=2,MRSE
      FV10 = 0.
      FV20 = 0.
      FV30 = 0.
      FV40 = 0.
      JSUM = 0
      DO J=SJB,SJE!1,MT-1
      IF (J.GE.MT) J=J-MT+1
      
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      ST= 0.5*FLOAT(MAX0(JP1-JM1,1))
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR= FLOAT(I-1)*DR

      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    VT2(I,J,K).NE.-999..and.VT2(IM1,J,K).NE.-999..and.
     &    VT2(IP1,J,K).NE.-999..and.VT2(I,JM1,K).NE.-999..and.
     &    VT2(I,JP1,K).NE.-999..and.VT2(I,J,KM1).NE.-999..and.
     &    VT2(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1

      ASFV= VT2(I,J,K)-BVR(I,K)
      ASFVIM1= VT2(IM1,J,K)-BVR(IM1,K)
      ASFVIP1= VT2(IP1,J,K)-BVR(IP1,K)
      ASFVKM1= VT2(I,J,KM1)-BVR(I,KM1)
      ASFVKP1= VT2(I,J,KP1)-BVR(I,KP1)
      DASFVDR= FNDIF1S(ASFVIM1,ASFV,ASFVIP1,DR1,DR2,SR)
      FV1(I,J,K)= -(VR2(I,J,K)-BUR(I,K))*
     &            (DASFVDR+(VT2(I,J,K)-BVR(I,K))/RR)
      FV10= FV10+ FV1(I,J,K)
      DASFVDZ= FNDIF1S(ASFVKM1,ASFV,ASFVKP1,DZ1,DZ2,SZ)
      FV2(I,J,K)= -(WRT(I,J,K)-BWR(I,K))* DASFVDZ
      FV20= FV20+ FV2(I,J,K)
      AP= PRT(I,J,K)-BPR(I,K)
      APJM1= PRT(I,JM1,K)-BPR(I,K)
      APJP1= PRT(I,JP1,K)-BPR(I,K)
      DAPDD= FNDIF1S(APJM1,AP,APJP1,DD1,DD2,ST)
      FV3(I,J,K)= (DRT(I,J,K)-BDR(I,K))*DAPDD/BDR(I,K)/BDR(I,K)
      FV30= FV30+ FV3(I,J,K)
      FV4(I,J,K)= RBLT2(I,J,K)
      FV40= FV40+ FV4(I,J,K)
      ENDIF
      ENDDO
      IF (JSUM.NE.0) THEN
      FV11(I,K)= FV10/FLOAT(JSUMALL(I,K)) !radial eddy angular momentum transport
      FV22(I,K)= FV20/FLOAT(JSUMALL(I,K)) !eddy vertical advection
      FV33(I,K)= FV30/FLOAT(JSUMALL(I,K)) !correlation between eddy pressure and density 
      FV44(I,K)= FV40/FLOAT(JSUMALL(I,K)) !turbulent momentum diffusion
CCC      FV11(I,K)= 0.
CCC      FV22(I,K)= 0.
CCC      FV33(I,K)= 0.
CCC      FV44(I,K)= 0.
      ELSE
      FV11(I,K)= 0.
      FV22(I,K)= 0.
      FV33(I,K)= 0.
      FV44(I,K)= 0.
      ENDIF
      ENDDO
      ENDDO

      DO K=1,NZ
      FV(1,K)=2*FV(2,K)-FV(3,K)
      FV11(1,K)=2*FV11(2,K)-FV11(3,K)
      FV22(1,K)=2*FV22(2,K)-FV22(3,K)
      FV33(1,K)=2*FV33(2,K)-FV33(3,K)
      FV44(1,K)=2*FV44(2,K)-FV44(3,K)
      ENDDO
      
c      DO K=16,NZ
c      DO I=1,MRSE
c      FV44(I,K)= 0.
c      FV22(I,K)= 0.
c      FV33(I,K)= 0.
c      ENDDO
c      ENDDO
    
      DO K=1,NZ
      DO I=1,MRSE
      FV(I,K)= FV11(I,K)+FV22(I,K)+FV33(I,K)+FV44(I,K)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      FVV(I,K)= CHI(I,K)*XII(I,K)*FV(I,K)
      ENDDO
      ENDDO
    
      DO K=1,NZ
      DO I=2,MRSE
c heat source
      HT10= 0.
      HT20= 0.
      HT30= 0.
      HT40= 0.
      HT50= 0.
      JSUM= 0
      DO J=SJB,SJE!1,MT-1
      IF (J.GE.MT) J=J-MT+1

      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      JM1= MAX0(J-1,1)
      JP1= MIN0(J+1,MT)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      ST= 0.5*FLOAT(MAX0(JP1-JM1,1))
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR= FLOAT(I-1)*DR

      IF (Z00(K).GE.terRT(I,J).and.Z00(K).GE.terRT(IM1,J).and.
     &    Z00(K).GE.terRT(IP1,J).and.Z00(K).GE.terRT(I,JM1).and.
     &    Z00(K).GE.terRT(I,JP1).and.Z00(KM1).GE.terRT(I,J).and.
     &    THETART(I,J,K).NE.-999..and.THETART(IM1,J,K).NE.-999..and.
     &    THETART(IP1,J,K).NE.-999..and.THETART(I,JM1,K).NE.-999..and.
     &    THETART(I,JP1,K).NE.-999..and.THETART(I,J,KM1).NE.-999..and.
     &    THETART(I,J,KP1).NE.-999.) THEN
      JSUM= JSUM+1

c asymmetric eddy heat source
      ASHT= THETART(I,J,K)-BTHETAR(I,K)
      ASHTIM1= THETART(IM1,J,K)-BTHETAR(IM1,K)
      ASHTIP1= THETART(IP1,J,K)-BTHETAR(IP1,K)
      ASHTJM1= THETART(I,JM1,K)-BTHETAR(I,K)
      ASHTJP1= THETART(I,JP1,K)-BTHETAR(I,K)
      ASHTKM1= THETART(I,J,KM1)-BTHETAR(I,KM1)
      ASHTKP1= THETART(I,J,KP1)-BTHETAR(I,KP1)
      DASHTDR= FNDIF1S(ASHTIM1,ASHT,ASHTIP1,DR1,DR2,SR)
      HT1(I,J,K)= -(VR2(I,J,K)-BUR(I,K))*DASHTDR
      HT10= HT10 + HT1(I,J,K)
      DASHTDD= FNDIF1S(ASHTJM1,ASHT,ASHTJP1,DD1,DD2,ST)
      HT2(I,J,K)= -(VT2(I,J,K)-BVR(I,K))*DASHTDD
      HT20= HT20 + HT2(I,J,K)
      DASHTDZ= FNDIF1S(ASHTKM1,ASHT,ASHTKP1,DZ1,DZ2,SZ)
      HT3(I,J,K)= -(WRT(I,J,K)-BWR(I,K))*DASHTDZ
      HT30= HT30 + HT3(I,J,K)
c      TT= THETART(I,J,K)/(P0/PRT(I,J,K))**0.286
      HT4(I,J,K)= HDIRT(I,J,K)
      HT40= HT40 + HT4(I,J,K)
      HT5(I,J,K)= RTHBLRT(I,J,K)
      HT50= HT50 + HT5(I,J,K)
      ENDIF
      ENDDO
      IF (JSUM.NE.0) THEN
      HT11(I,K) = HT10/FLOAT(JSUMALL(I,K)) !radial eddy temperature advection
      HT22(I,K) = HT20/FLOAT(JSUMALL(I,K)) !tangential eddy temperature advection
      HT33(I,K) = HT30/FLOAT(JSUMALL(I,K)) !vertical eddy temperature advection
      HT44(I,K) = HT40/FLOAT(JSUMALL(I,K)) !diabatic heating
      HT55(I,K) = HT50/FLOAT(JSUMALL(I,K)) !turbulent heat diffusion
CCC      HT11(I,K)= 0.
CCC      HT22(I,K)= 0.
CCC      HT33(I,K)= 0.
CCC      HT44(I,K)= 0.
CCC      HT55(I,K)= 0.
      ELSE
      HT11(I,K)= 0.
      HT22(I,K)= 0.
      HT33(I,K)= 0.
      HT44(I,K)= 0.
      HT55(I,K)= 0.
      ENDIF
      ENDDO
      ENDDO
      
c      DO K=1,75
c      DO I=1,MRSE
c      HT11(I,K)= 0.
c      HT22(I,K)= 0.
c      HT33(I,K)= 0.
c      ENDDO
c      ENDDO   
c      Print*,'HT3',HT11(3,3)
c      Print*,'HT65',HT11(65,65)
      DO K=1,NZ
      HT(1,K)=2*HT(2,K)-HT(3,K)
      HT11(1,K)=2*HT11(2,K)-HT11(3,K)
      HT22(1,K)=2*HT22(2,K)-HT22(3,K)
      HT33(1,K)=2*HT33(2,K)-HT33(3,K)
      HT44(1,K)=2*HT44(2,K)-HT44(3,K)
      HT55(1,K)=2*HT55(2,K)-HT55(3,K)
      ENDDO
      DO K=1,NZ
      DO I=1,MRSE
      HT(I,K) = HT11(I,K)+HT22(I,K)+HT33(I,K)+HT44(I,K)+HT55(I,K)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=1,MRSE
      HTDR(I,K)=GG*CHI(I,K)*CHI(I,K)*(1.+epsi*BQVR(I,K))*HT(I,K)
      HTDZ(I,K)=CC(I,K)*CHI(I,K)*CHI(I,K)*(1.+epsi*BQVR(I,K))*HT(I,K)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR= FLOAT(I-1)*DR
      FVtal(I,K)= -FNDIF1S(FVV(I,KM1),FVV(I,K),FVV(I,KP1),DZ1,DZ2,SZ)
      DHTDR(I,K)= FNDIF1S(HTDR(IM1,K),HTDR(I,K),HTDR(IP1,K),DR1,DR2,SR)
      DHTDZ(I,K)= FNDIF1S(HTDZ(I,KM1),HTDZ(I,K),HTDZ(I,KP1),DZ1,DZ2,SZ)
      ENDDO
      ENDDO

      DO K=1,NZ
      FVtal(1,K)= FVtal(2,K)
      DHTDR(1,K)= DHTDR(2,K)
      DHTDZ(1,K)= DHTDZ(2,K)
      ENDDO
      
c-----------------------------------------------------------------
c calculate Vdot
      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      DVDT= (AVTSE(I,K)-AVTSEb(I,K))/DT
      VDO2= BUR(I,K)*ETA(I,K)
      DVDZ= FNDIF1S(AVTSE(I,KM1),AVTSE(I,K),AVTSE(I,KP1),
     $  DZ1,DZ2,SZ)
      VDO3 = BWR(I,K)*DVDZ
      VDOT(I,K)= DVDT+VDO2+VDO3-FV(I,K)
      VDDZ(I,K)= CHI(I,K)*XII(I,K)*VDOT(I,K)
      ENDDO
      ENDDO

      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR= FLOAT(I-1)*DR
      DVDDZ(I,K)= FNDIF1S(VDDZ(I,KM1),VDDZ(I,K),VDDZ(I,KP1),
     $  DZ1,DZ2,SZ)
CCC      DVDDZ(I,K)= 0.
      ENDDO
      ENDDO
      DO K=1,NZ
      DVDDZ(1,K)=2*DVDDZ(2,K)-DVDDZ(3,K)
      ENDDO
c--------------------------------------------------------------

      DO K=1,NZ
      DO I=1,MRSE
c      FVtal(I,K)= 0.
c      DHTDR(I,K)= 0.
c      DHTDZ(I,K)= 0.
      SS(I,K)= FVtal(I,K)+DHTDR(I,K)+DHTDZ(I,K)+
     $  DUDDT(I,K)-DWDDT(I,K)-DVDDZ(I,K)
      SS(I,K)= SS(I,K)
      ENDDO
      ENDDO
      call MAXMIN2("SS",SS,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("FVtal",FVtal,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("DHTDR",DHTDR,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("DHTDZ",DHTDZ,MRSE,1,NZ,MR,NZ,0,0)
      
      DO K=1,NZ
      DO I=1,MRSE
      WK2D(I,K,1)= VPT(I,K)
      WK2D(I,K,2)= AVTSE(I,K)
      WK2D(I,K,3)= ATI(I,K)
      WK2D(I,K,4)= BTI1(I,K)
      WK2D(I,K,5)= CTI(I,K)
      WK2D(I,K,6)= FV(I,K)
      WK2D(I,K,7)= FVV(I,K)
      WK2D(I,K,8)= HT(I,K)
      WK2D(I,K,9)= HTDR(I,K)
      WK2D(I,K,10)= HTDZ(I,K)
      WK2D(I,K,11)= CON2(I,K)
      WK2D(I,K,12)= BTI2(I,K)
      WK2D(I,K,13)= PII(I,K)
      WK2D(I,K,14)= CON1(I,K)
      WK2D(I,K,15)= CHI(I,K)
      WK2D(I,K,16)= CC(I,K)
      WK2D(I,K,17)= CHIC(I,K)
      WK2D(I,K,18)= XII(I,K)
      WK2D(I,K,19)= ETA(I,K)
      WK2D(I,K,20)= HT11(I,K)
      WK2D(I,K,21)= HT22(I,K)
      WK2D(I,K,22)= HT33(I,K)
      WK2D(I,K,23)= HT44(I,K)
      WK2D(I,K,24)= HT55(I,K)
      WK2D(I,K,25)= FV11(I,K)
      WK2D(I,K,26)= FV22(I,K)
      WK2D(I,K,27)= FV33(I,K)
      WK2D(I,K,28)= FV44(I,K)
      WK2D(I,K,29)= FVtal(I,K)
      WK2D(I,K,30)= DHTDR(I,K)
      WK2D(I,K,31)= DHTDZ(I,K)
      WK2D(I,K,32)= VDOT(I,K)
      WK2D(I,K,33)= SS(I,K)
      WK2D(I,K,34)= PVRT1(I,K)
      WK2D(I,K,35)= PVRT2(I,K)
      WK2D(I,K,36)= PVRT(I,K)
      WK2D(I,K,37)= NRO(I,K)
      WK2D(I,K,38)= DVDDZ(I,K)
      WK2D(I,K,39)= WDOT(I,K)
      WK2D(I,K,40)= UDOT(I,K)
      WK2D(I,K,41)= DUDDT(I,K)
      WK2D(I,K,42)= DWDDT(I,K)
      ENDDO
      ENDDO
      print*,'VPT44',VPT(50,50)

      open(unit=44,file=trim(filename(FN))//'_abc.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MR*NZ*4)
      DO L=1,33
      WRITE(44,rec=L) WK2D(:,:,L)
      ENDDO
      WRITE(44,rec=34) WK2D(:,:,38)
      WRITE(44,rec=35) WK2D(:,:,39)
      WRITE(44,rec=36) WK2D(:,:,40)
      WRITE(44,rec=37) WK2D(:,:,41)
      WRITE(44,rec=38) WK2D(:,:,42)
      CLOSE (UNIT=44)
      
      open(unit=45,file=trim(filename(FN))//'_PV.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MR*NZ*4)
      DO L=34,37
      WRITE(45,rec=L-33) WK2D(:,:,L)
      ENDDO
      CLOSE (UNIT=45)

      open(unit=46,file=trim(filename(FN))//'_3D.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MR*MT*NZ*4)
      WRITE(46,rec=1) FV1(:,:,:)
      WRITE(46,rec=2) FV2(:,:,:)
      WRITE(46,rec=3) FV3(:,:,:)
      WRITE(46,rec=4) FV4(:,:,:)
      WRITE(46,rec=5) HT1(:,:,:)
      WRITE(46,rec=6) HT2(:,:,:)
      WRITE(46,rec=7) HT3(:,:,:)
      WRITE(46,rec=8) HT4(:,:,:)
      WRITE(46,rec=9) HT5(:,:,:)
      CLOSE (UNIT=46)
c      
      ENDIF
c-------------------------------------------------------
! begin integration
       IF (INTE.eq.1) THEN
       open(unit=55,file=trim(filename(FN))//'_abc.dat',
     &   form='unformatted',
     &   status='unknown',access='stream')
         read(55) ((VPT(i,k),i=1,MR),k=1,NZ)
         read(55) ((AVTSE(i,k),i=1,MR),k=1,NZ)
         read(55) ((ATI(i,k),i=1,MR),k=1,NZ)
         read(55) ((BTI1(i,k),i=1,MR),k=1,NZ)
         read(55) ((CTI(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV(i,k),i=1,MR),k=1,NZ)
         read(55) ((FVV(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT(i,k),i=1,MR),k=1,NZ)
         read(55) ((HTDR(i,k),i=1,MR),k=1,NZ)
         read(55) ((HTDZ(i,k),i=1,MR),k=1,NZ)
         read(55) ((CON2(i,k),i=1,MR),k=1,NZ)
         read(55) ((BTI2(i,k),i=1,MR),k=1,NZ)
         read(55) ((PII(i,k),i=1,MR),k=1,NZ)
         read(55) ((CON1(i,k),i=1,MR),k=1,NZ)
         read(55) ((CHI(i,k),i=1,MR),k=1,NZ)
         read(55) ((CC(i,k),i=1,MR),k=1,NZ)
         read(55) ((CHIC(i,k),i=1,MR),k=1,NZ)
         read(55) ((XII(i,k),i=1,MR),k=1,NZ)
         read(55) ((ETA(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT11(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT22(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT33(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT44(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT55(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV11(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV22(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV33(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV44(i,k),i=1,MR),k=1,NZ)
         read(55) ((FVtal(i,k),i=1,MR),k=1,NZ)
         read(55) ((DHTDR(i,k),i=1,MR),k=1,NZ)
         read(55) ((DHTDZ(i,k),i=1,MR),k=1,NZ)
         read(55) ((VDOT(i,k),i=1,MR),k=1,NZ)
         read(55) ((SS(i,k),i=1,MR),k=1,NZ)
       CLOSE (UNIT=55)

      print*,'ATI',ATI(3,3)
      Print*,'AVTSE',AVTSE(3,3)
      DO K=1,NZ
      DO I=1,MRSE
      DSEST1(I,K)=0.
      DSEST2(I,K)=0.
      ENDDO
      ENDDO

      WSO= 0.8 
      NTOTAL=10000000
      NSEE= 1
      EPS=1.E-20
      DRD= DR
      WSOD= WSO
c
      DO 57 NN= 1,NTOTAL
c add SOR
      IF (IDZCONST .EQ. 1) THEN
      DZD= DZ(1)
      DO K=2,NZ-1
      DO I=2,MRSE-1
      DGAM(I,K)= 2.*ATI(I,K)/(DRD*DRD)+2.*CTI(I,K)/(DZD*DZD)
      ENDDO
      ENDDO
      DO K=2,NZ-1
      DO I=2,MRSE-1
      AIP12= (ATI(I+1,K)+ATI(I,K))/2.
      AIM12= (ATI(I,K)+ATI(I-1,K))/2.
      CKP12= (CTI(I,K+1)+CTI(I,K))/2.
      CKM12= (CTI(I,K)+CTI(I,K-1))/2.
      BTI1D1= BTI1(I+1,K)
      BTI1D2= BTI1(I-1,K)
      IF (I.EQ.2) THEN
      AIM12= A1P(K)
      ENDIF
      DSS= SS(I,K)
      DA11(I,K)= AIP12*(DSEST1(I+1,K)-DSEST1(I,K))/DRD
      DA12(I,K)= AIM12*(DSEST1(I,K)-DSEST1(I-1,K))/DRD
      DB11(I,K)= BTI1D1*(DSEST1(I+1,K+1)-DSEST1(I+1,K-1))/(2.*DZD)
      DB12(I,K)= BTI1D2*(DSEST1(I-1,K+1)-DSEST1(I-1,K-1))/(2.*DZD)
      DC11(I,K)= CKP12*(DSEST1(I,K+1)-DSEST1(I,K))/DZD
      DC12(I,K)= CKM12*(DSEST1(I,K)-DSEST1(I,K-1))/DZD
      DB21(I,K)= BTI2(I,K+1)*(DSEST1(I+1,K+1)-DSEST1(I-1,K+1))/(2.*DR)
      DB22(I,K)= BTI2(I,K-1)*(DSEST1(I+1,K-1)-DSEST1(I-1,K-1))/(2.*DR)
      DSERSO1(I,K)= (DA11(I,K)-DA12(I,K))/DRD
      DSERSO2(I,K)= (DB11(I,K)-DB12(I,K))/(2.*DRD)
      DSERSO3(I,K)= (DC11(I,K)-DC12(I,K))/DZD
      DSERSO4(I,K)= (DB21(I,K)-DB22(I,K))/(2.*DZD)
      DSERSO(I,K)= DSERSO1(I,K)+DSERSO2(I,K)+DSERSO3(I,K)+DSERSO4(I,K)
     $  -DSS
      DSEST2(I,K)= DSEST1(I,K)+WSOD*(1./DGAM(I,K))*DSERSO(I,K)
      SERSO(I,K)= DSERSO(I,K)
      ENDDO
      ENDDO
      ELSE
      DO K=1,NZ
      DZDP(K)= DZ(K)
      ENDDO
      DO K=2,NZ-1
      DO I=2,MRSE-1
      CKP12= (CTI(I,K+1)+CTI(I,K))/2.
      CKM12= (CTI(I,K)+CTI(I,K-1))/2.
      DGAM(I,K)=(2.*ATI(I,K)/(DRD*DRD)+
     $  CKP12/(DZDP(K)*(DZDP(K)+DZDP(K-1))/2.)+
     $  CKM12/(DZDP(K-1)*(DZDP(K)+DZDP(K-1))/2.))
      ENDDO
      ENDDO
      DO K=2,NZ-1
      DO I=2,MRSE-1
      AIP12= (ATI(I+1,K)+ATI(I,K))/2.
      AIM12= (ATI(I,K)+ATI(I-1,K))/2.
      CKP12= (CTI(I,K+1)+CTI(I,K))/2.
      CKM12= (CTI(I,K)+CTI(I,K-1))/2.
      BTI1D1= BTI1(I+1,K)
      BTI1D2= BTI1(I-1,K)
      IF (I.EQ.2) THEN
      AIM12= A1P(K)
      ENDIF
      DSS= SS(I,K)
      DA11(I,K)= AIP12*(DSEST1(I+1,K)-DSEST1(I,K))/DRD
      DA12(I,K)= AIM12*(DSEST1(I,K)-DSEST1(I-1,K))/DRD
      DB11(I,K)= BTI1D1*(DSEST1(I+1,K+1)-DSEST1(I+1,K-1))
     $  /(DZDP(K)+DZDP(K-1))
      DB12(I,K)= BTI1D2*(DSEST1(I-1,K+1)-DSEST1(I-1,K-1))
     $  /(DZDP(K)+DZDP(K-1))
      DC11(I,K)= CKP12*(DSEST1(I,K+1)-DSEST1(I,K))/DZDP(K)
      DC12(I,K)= CKM12*(DSEST1(I,K)-DSEST1(I,K-1))/DZDP(K-1)
      DB21(I,K)= BTI2(I,K+1)*(DSEST1(I+1,K+1)-DSEST1(I-1,K+1))/(2.*DR)
      DB22(I,K)= BTI2(I,K-1)*(DSEST1(I+1,K-1)-DSEST1(I-1,K-1))/(2.*DR)
      DSERSO1(I,K)= (DA11(I,K)-DA12(I,K))/DRD
      DSERSO2(I,K)= (DB11(I,K)-DB12(I,K))/(2.*DRD)
      DSERSO3(I,K)= (DC11(I,K)-DC12(I,K))/((DZDP(K)+DZDP(K-1))/2.)
      DSERSO4(I,K)= (DB21(I,K)-DB22(I,K))/(DZDP(K)+DZDP(K-1))
      DSERSO(I,K)= DSERSO1(I,K)+DSERSO2(I,K)+DSERSO3(I,K)+DSERSO4(I,K)
     $  -DSS
      DSEST2(I,K)= DSEST1(I,K)+WSOD*(1./DGAM(I,K))*DSERSO(I,K)
      SERSO(I,K)= DSERSO(I,K)
      ENDDO
      ENDDO
      ENDIF
c
      SERSOmax=0.
      DO K=2,NZ-1
      DO I=2,MRSE-1
      SERSOmax= AMAX1(SERSOmax,abs(SERSO(I,K)))
      ENDDO
      ENDDO
c      call MAXMIN2("GAM",DGAM,MRSE,2,NZ-1,MR,NZ,1,0)
c      call MAXMIN2("SERSO",SERSO,MRSE,2,NZ-1,MR,NZ,1,0)
c      call MAXMIN2("DSERSO1",DSERSO1,MRSE,2,NZ-1,MR,NZ,1,0)
c      call MAXMIN2("DSERSO2",DSERSO2,MRSE,2,NZ-1,MR,NZ,1,0)
c      call MAXMIN2("DSERSO3",DSERSO3,MRSE,2,NZ-1,MR,NZ,1,0)
c      call MAXMIN2("DSERSO4",DSERSO4,MRSE,2,NZ-1,MR,NZ,1,0)
! set boundary conditions
      DO K=1,NZ
      DSEST2(1,K)= DSEST2(2,K)
      DSEST2(1,K)= 0.
      DSEST2(MRSE,K)= DSEST2(MRSE-1,K)
      ENDDO         
      DO I=1,MRSE
      DSEST2(I,1)= 0.
      DSEST2(I,NZ)= DSEST2(I,NZ-1)
      DSEST2(I,NZ)= 0.
      ENDDO
      DO K=1,NZ
      DO I=1,MRSE      
      SEST1(I,K)= DSEST1(I,K)
      SEST2(I,K)= DSEST2(I,K)
      ENDDO
      ENDDO
      SEST2max=0.
      DO K=2,NZ-1
      DO I=2,MRSE-1
      SEST2max= AMAX1(SEST2max,abs(SEST2(I,K)))
      ENDDO
      ENDDO
      IF (mod(NN,NSEE).EQ.0) THEN 
      PRINT *,'NN=',NN,'MAX |SERSO|=',SERSOmax
      PRINT *,'SEST2max',SEST2max
      call MAXMIN2("SEST1",SEST1,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("SEST2",SEST2,MRSE,1,NZ,MR,NZ,0,0)
      call MAXMIN2("SERSO",SERSO,MRSE,2,NZ-1,MR,NZ,1,0)
      ENDIF
      DO K=1,NZ
      DO I=1,MRSE
      SEST1(I,K)= SEST2(I,K)
      DSEST1(I,K)= DSEST2(I,K)
      ENDDO
      ENDDO

      open(unit=59,file=trim(pgfile)//'_cyli2_1.dat',form=
     & 'unformatted',status='unknown',access='stream')
         read(59) ((BUR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BVR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BWR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BPR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BDR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BHDIR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BTHETAR(i,k),i=1,MR),k=1,NZ)
         read(59) ((BQVR(i,k),i=1,MR),k=1,NZ)
         read(59) ((FRT(i,k),i=1,MR),k=1,NZ)
      CLOSE (UNIT=59)
      DO K=1,NZ
      DO I=2,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR= FLOAT(I-1)*DR
c transverse flow
      DSEDZ(I,K)= FNDIF1S(SEST2(I,KM1),SEST2(I,K),SEST2(I,KP1),
     & DZ1,DZ2,SZ) 
      DSEDR(I,K)= FNDIF1S(SEST2(IM1,K),SEST2(I,K),SEST2(IP1,K),
     & DR1,DR2,SR)
      SEUU(I,K)= -DSEDZ(I,K)/(RR*BDR(I,K))
      SEWW(I,K)= DSEDR(I,K)/(RR*BDR(I,K))
      ENDDO
      ENDDO
      DO K=1,NZ
      SEUU(1,K)= 0.
      SEWW(1,K)= 0.
      ENDDO
      UUmax=0.
      WWmax=0.
      DO K=1,NZ
      DO I=2,MRSE-1
      UUmax= AMAX1(UUmax,abs(SEUU(I,K)))
      WWmax= AMAX1(WWmax,abs(SEWW(I,K)))
      ENDDO
      ENDDO
      IF (mod(NN,NSEE).EQ.0) THEN
      PRINT*,'UUmax',UUmax,'WWmax',WWmax
      call MAXMIN2("SEUU",SEUU,MRSE,1,NZ,MR,NZ,0,1)
      call MAXMIN2("SEWW",SEWW,MRSE,1,NZ,MR,NZ,0,1)
      call MAXMIN2("DSEDZ",DSEDZ,MRSE,1,NZ,MR,NZ,0,1)
      call MAXMIN2("DSEDR",DSEDR,MRSE,1,NZ,MR,NZ,0,1)
      PRINT*,'SEUU(103,3)',SEUU(103,3)
      PRINT*,'DSEDZ(103,3)',DSEDZ(103,3)
      ENDIF

      IF (SERSOmax.LE.EPS .OR. SERSOmax  .GE. 1.E10) go to 58
            
 57   CONTINUE
 58   CONTINUE
c      DO K=1,NZ
c      DO I=1,MRSE
c      WK2DSE(I,K,1)= WK2DSE(I,K,1) + SEUU(I,K)
c      WK2DSE(I,K,2)= WK2DSE(I,K,2) + SEWW(I,K)
c      ENDDO
c      ENDDO

c      open(unit=66,file=trim(pgfile)//'_SE.dat',form='unformatted',
c     &  access='direct',status='unknown',recl=MRSE*NZ*4)
c      DO L=1,3
c      WRITE(66,rec=L) WK2DSE(:,:,L)
c      ENDDO
c      CLOSE (UNIT=66)
c  azimuthal mean tangential wind
      open(unit=55,file=trim(filename(FN))//'_abc.dat',
     &  form='unformatted',
     &  status='unknown',access='stream')
         read(55) ((VPT(i,k),i=1,MR),k=1,NZ)
         read(55) ((AVTSE(i,k),i=1,MR),k=1,NZ)
         read(55) ((ATI(i,k),i=1,MR),k=1,NZ)
         read(55) ((BTI1(i,k),i=1,MR),k=1,NZ)
         read(55) ((CTI(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV(i,k),i=1,MR),k=1,NZ)
         read(55) ((FVV(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT(i,k),i=1,MR),k=1,NZ)
         read(55) ((HTDR(i,k),i=1,MR),k=1,NZ)
         read(55) ((HTDZ(i,k),i=1,MR),k=1,NZ)
         read(55) ((CON2(i,k),i=1,MR),k=1,NZ)
         read(55) ((BTI2(i,k),i=1,MR),k=1,NZ)
         read(55) ((PII(i,k),i=1,MR),k=1,NZ)
         read(55) ((CON1(i,k),i=1,MR),k=1,NZ)
         read(55) ((CHI(i,k),i=1,MR),k=1,NZ)
         read(55) ((CC(i,k),i=1,MR),k=1,NZ)
         read(55) ((CHIC(i,k),i=1,MR),k=1,NZ)
         read(55) ((XII(i,k),i=1,MR),k=1,NZ)
         read(55) ((ETA(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT11(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT22(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT33(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT44(i,k),i=1,MR),k=1,NZ)
         read(55) ((HT55(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV11(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV22(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV33(i,k),i=1,MR),k=1,NZ)
         read(55) ((FV44(i,k),i=1,MR),k=1,NZ)
         read(55) ((FVtal(i,k),i=1,MR),k=1,NZ)
         read(55) ((DHTDR(i,k),i=1,MR),k=1,NZ)
         read(55) ((DHTDZ(i,k),i=1,MR),k=1,NZ)
         read(55) ((VDOT(i,k),i=1,MR),k=1,NZ)
         read(55) ((SS(i,k),i=1,MR),k=1,NZ)
      close(UNIT=55)
      Print*,'AVTSE',AVTSE(3,3)

      DO K=1,NZ
      DO I=1,MRSE
      IM1= MAX0(I-1,1)
      IP1= MIN0(I+1,MRSE)
      KM1= MAX0(K-1,1)
      KP1= MIN0(K+1,NZ)
      SR= 0.5*FLOAT(MAX0(IP1-IM1,1))
      SZ= 0.5*FLOAT(MAX0(KP1-KM1,1))
      DR1= DR
      DR2= DR
      DD1= DD(I)
      DD2= DD(I)
      DZ1= DZ(KM1)
      DZ2= DZ(K)
      RR = FLOAT(I-1)*DR
c gradient-wind balance
      AGV1(I,K) = -SEUU(I,K)*ETA(I,K)
      DAGVDZ =  FNDIF1S(AVTSE(I,KM1),AVTSE(I,K),AVTSE(I,KP1),
     $  DZ1,DZ2,SZ)
      AGV2(I,K) = -SEWW(I,K)*DAGVDZ
      TWT(I,K) = FV(I,K)+AGV1(I,K)+AGV2(I,K)
c      TWT(I,K) = +AGV1+AGV2
      ENDDO
      ENDDO
      
c      WK2D(I,K,24)= WK2D(I,K,24) + TWT(I,K)

      DO K=1,NZ
      DO I=1,MRSE
      WK2DSE(I,K,1)= SEUU(I,K)
      WK2DSE(I,K,2)= SEWW(I,K)
      WK2DSE(I,K,3)= TWT(I,K)
      WK2DSE(I,K,4)= SEST2(I,K)
      WK2DSE(I,K,5)= FV(I,K)
      WK2DSE(I,K,6)= AGV1(I,K)
      WK2DSE(I,K,7)= AGV2(I,K)
      ENDDO
      ENDDO
      call MAXMIN2("TWT",TWT,MRSE,1,NZ,MR,NZ,0,0)
      open(unit=66,file=trim(filename(FN))//'_SE.dat',
     &  form='unformatted',
     &  access='direct',status='unknown',recl=MRSE*NZ*4)
      DO L=1,7
      WRITE(66,rec=L) WK2DSE(:,:,L)
      ENDDO
      CLOSE (UNIT=66)
      ENDIF
      IF (ISHEAR.EQ.1) GO TO 302
  303 CONTINUE
c      GO TO 88
c 99   CONTINUE
c
      DEALLOCATE(theta,pre,den,uu,vv,
     &  ww,qv,qc,qr,qi,qs,qg,rubl,rvbl,rucu,rvcu,rthbl,
     &  rthcu,rthlw,rthsw,rthdia,
     &  urt,vrt,wrt,prt,drt,
     &  vr2,vt2,
     &  HDIRT,THETART,QVRT,
     &  BHDIR,BTHETAR,BQVR,
     &  PII,VPT,CHI,AVTSE,ATI,CC,CHIC,BTI1,BTI2,XII,ETA,CTI,CON1,CON2,
     &  FV,HT,DSEST,SEST1,SEST2,SEDR,SEDZ,FVV,HTDR,HTDZ,
     &  SEUU,SEWW,DSEDZ,DSEDR,RBDR)

      STOP
      END
C
C--------------------------------------------------------------------
       FUNCTION FNDIF1S(Q1,Q2,Q3,D1,D2,SZ)
C
       DD1= D1*D1
       DD2= D2*D2
       QDX= (DD1*Q3-DD2*Q1-DD1*Q2+DD2*Q2)/(SZ*(DD1*D2+DD2*D1))
C
       FNDIF1S= QDX
C
       RETURN
       END

c--------------------------------------------------------------------
       FUNCTION FNDIF2(QM,Q,QP,DM,D)
C
       DD= 0.5*(DM+D)
       QDXX= ( (QP-Q)/D - (Q-QM)/DM )/DD

       FNDIF2= QDXX
C
       RETURN
       END

c--------------------------------------------------------------------
       FUNCTION FNDIFXY(QMM,QPM,QMP,QPP,D1M,D1,D2M,D2)
C
       DD= (D1+D1M)*(D2+D2M)
       QDXY= (QPP+QMM-QMP-QPM)/DD

       FNDIFXY= QDXY
C
       RETURN
       END

C---------------------------------------------------------------------------------
      SUBROUTINE VINP(WW,XX,WW2,NX,NY,NZ,NZ2,W1D,X1D,W1D2,X2)
C
c WW: 3D values at XX grids (NX,NY,NZ)
c WW2: interpolated values on new vertical grids X2(NZ2)

      REAL*8  WW(NX,NY,NZ),XX(NX,NY,NZ)
      REAL    WW2(NX,NY,NZ2),W1D(NZ),X1D(NZ)
      REAL    W1D2(NZ2),X2(NZ2)
C
      do i=1,NX
      do j=1,NY
      do k=1,NZ
      w1d(k)= WW(i,j,k)
      x1d(k)= XX(i,j,k)
      enddo
cc      print *,'w1d=',w1d
cc      print *,'x1d=',x1d
      do k=1,NZ2
      xp= X2(k)
      w1d2(k)= CUBLA(xp,NZ,x1d,w1d,CHECK) 
      enddo
cc      print *,'w1d2=',w1d2
      do k=1,NZ2
      WW2(i,j,k)= w1d2(k)
      enddo
      enddo
      enddo

      return 
      end 
C
C------------------------------------------------------------------
       FUNCTION CUBLA(XP,NX,TDX,Q,CHECK)
       DIMENSION TDX(NX),Q(NX)
C
C  for stretched grids with unknown XP
c
       CHECK= 0.
C
       N= NX
       X= XP
C
C  iIVC= 1: upward increase; IVC=2: upward decrease
C  check coordinate
       IF (TDX(NX).GT.TDX(1)) then
         IVC= 1 
       ELSE
         IVC= 2
       ENDIF
C
C  determine I for X
C
       IF (IVC.EQ.1) THEN
       IF (X.LE.TDX(1)) THEN
         IIP= 1
         GO TO 10
       ENDIF
       IF (X.GE.TDX(N)) THEN
         IIP= N-1
         GO TO 10
       ENDIF
       ELSE
       IF (X.GE.TDX(1)) THEN
         IIP= 1
         GO TO 10
       ENDIF
       IF (X.LE.TDX(N)) THEN
         IIP= N-1
         GO TO 10
       ENDIF
       ENDIF

       IF (IVC.EQ.1) THEN
       DO II=1,N-1
       IF (X.GE.TDX(II).AND.X.LT.TDX(II+1)) THEN
         IIP= II
         GO TO 10
       ENDIF
       ENDDO
       ELSE
       DO II=1,N-1
       IF (X.LE.TDX(II).AND.X.GT.TDX(II+1)) THEN
         IIP= II
         GO TO 10
       ENDIF
       ENDDO
       ENDIF
C
  10   IIP= MIN0(IIP,N)
C
       I= IIP
C
       IF (I.GE.2.AND.I.LE.N-2) THEN
C
       IM1= MAX0(I-1,1)
       IP1= MIN0(I+1,N)
       IP2= MIN0(I+2,N)
       X0= TDX(IM1)
       X1= TDX(I)
       X2= TDX(IP1)
       X3= TDX(IP2)
C
       XX1= X - X1
       XX2= X - X2
       XX3= X - X3
       XX0= X - X0
       X10= X1 - X0
       X20= X2 - X0
       X30= X3 - X0
       X21= X2 - X1
       X32= X3 - X2
       X31= X3 - X1
C
       P3X= -(XX1/X10)*(XX2/X20)*(XX3/X30)*Q(IM1)
     $   +  (XX0/X10)*(XX2/X21)*(XX3/X31)*Q(I)
     $   - (XX0/X20)*(XX1/X21)*(XX3/X32)*Q(IP1)
     $   + (XX0/X30)*(XX1/X31)*(XX2/X32)*Q(IP2)
C
       ELSE

       IP1= MAX0(I+1,1)
       P3X= ( (X-TDX(I))*Q(IP1) + (TDX(IP1)-X)*Q(I) )/(TDX(IP1)-TDX(I))

       ENDIF
C
       CUBLA= P3X
C
       RETURN
       END

C
C--------------------------------------------------------------------------
       SUBROUTINE ASYMV(UXY,VXY,NX,NY,TDX,TDY,XCS,YCS,
     &   IMMCH,SPV,IREMOVE,terRT0,terrain0)
       PARAMETER (MR=400,MT=361,MRF=6)  !'CHANGE ME'
       DOUBLE PRECISION XCS,YCS,DIS,XXR,X1,Y1,XC,YC
       DIMENSION UXY(NX,NY),VXY(NX,NY),TDX(NX),TDY(NY)
       DIMENSION URT(MR,MT),VRT(MR,MT),VR2(MR,MT),VT2(MR,MT)
       DIMENSION BUV(MR),BVR(MR),BVT(MR),XXR(MRF),FFR(MRF)
       DIMENSION terRT0(MR,MT),terrain0(NX,NY)
       DATA PI/3.14159265358979/
C
C  finds the symmetric center of a wind field  and removes the symmetric
C  part
C
C  ON RETURN: 
C      UXY,VXY:  with the symmetric component removed and filled beyond
C      MR
C                filled with SPV beyond MR
C  TDX, TDY: X,Y coordinates (monotonically increase)
C  SPV= 0 to set wind vectors beyond MR
C  XCS,YCS:  > 0 used as the coordinate for neibourhood searching; 
C            < 0 used as the grid index for neibourhood searching for
C            the center
C            -9999 for global seraching of the center depending on IMMCH
C  IMMCH= 1 (0) searching maximum (minimum) defined as the center
C  Note:  Set MR=NX and MT=NY or so
C
       XC= XCS
       YC= YCS
       DR= ABS(TDX(2)-TDX(1))
       ICHECK= 0
C
       CALL CYLIN(UXY,NX,NY,TDX,TDY,URT,MR,MT,BUV,DR,XC,YC,SPV,IMMCH,
     &            0,0,terRT0(:,:),0.,0.,0,terrain0(:,:),0,terrain0(:,:),
     &            terrain0(:,:),0)
       CALL CYLIN(VXY,NX,NY,TDX,TDY,VRT,MR,MT,BUV,DR,XC,YC,SPV,IMMCH,
     &            0,0,terRT0(:,:),0.,0.,0,terrain0(:,:),0,terrain0(:,:),
     &            terrain0(:,:),0)
C
       CALL SVRVT(URT,VRT,VR2,VT2,MR,MT,1)
C
       XCS=XC
       YCS=YC
C
C
c  the minimum radius for use
       TRX1= ABS(TDX(1)-XC)
       TRX2= ABS(TDX(NX)-XC)
       TRY1= ABS(TDY(1)-YC)
       TRY2= ABS(TDY(NY)-YC)
       TRMN= AMIN1(TRX1,TRX2,TRY1,TRY2)
       IRMN= INT(TRMN/DR)+1
c       print *,'Minimum radius dimension IRMN= ',IRMN
       IF (IRMN.GE.MR) THEN
c       print *,'!!  Need to increase MR...'
       ENDIF
C
       IF (ICHECK.GE.2) THEN
       DO J=1,1
c       print *,'J=',J,'  VR2=',(VR2(I,J),I=1,IRMN)
       ENDDO
       DO J=1,1
c       print *,'J=',J,'  VT2=',(VT2(I,J),I=1,IRMN)
       ENDDO
       ENDIF
C
       DO I=1,IRMN
       TVR2= 0.
       TVT2= 0.
       JSUM1= 0
       JSUM2= 0
       DO J=1,MT-1
       IF (VR2(I,J).NE.SPV) THEN
       JSUM1= JSUM1+1
       TVR2= TVR2+ VR2(I,J)
       ENDIF
       IF (VT2(I,J).NE.SPV) THEN
       JSUM2= JSUM2+1
       TVT2= TVT2+ VT2(I,J)
       ENDIF
       ENDDO
       BVR(I)= TVR2/FLOAT(JSUM1)
       BVT(I)= TVT2/FLOAT(JSUM2)
       ENDDO
       IF (ICHECK.GE.1) THEN
c       print *,'BVR=',(BVR(I),I=1,IRMN)
c       print *,'BVT=',(BVT(I),I=1,IRMN)
       ENDIF
C
C  Use the FT method to get wavenumber-N components
       ISYM= 1
cc       IF (ISYM.EQ.0) THEN
cc       DO J=1,MT
cc       DO I=1,MR
cc       VR2(I,J)= VR2(I,J)-BVR(I)
cc       VT2(I,J)= VT2(I,J)-BVT(I)
cc       ENDDO
cc       ENDDO
cc       ENDIF
       IF (ISYM.EQ.1.AND.IREMOVE.NE.0) THEN
       DO I=1,MR
       CALL FT(VR2(I,1:MT-1),1,1,1,MT-1) 
       CALL FT(VT2(I,1:MT-1),1,1,1,MT-1) 
       ENDDO
       DO I=1,MR
       VR2(I,MT)= VR2(I,1)
       VT2(I,MT)= VT2(I,1)
       ENDDO
C  UXY used to record VRXY, VXY used to record VTXY
       CALL RTTOXY(UXY,NX,NY,TDX,TDY,VR2,MR,MT,BUV,DR,XC,YC,SPV,IMMCH)
       CALL RTTOXY(VXY,NX,NY,TDX,TDY,VT2,MR,MT,BUV,DR,XC,YC,SPV,IMMCH)
       DO J=1,NY
       DO I=1,NX
       X1= TDX(I)-XC
       Y1= TDY(J)-YC
       DIS= SQRT(X1*X1+Y1*Y1)
       VRX= UXY(I,J)
       VTX= VXY(I,J)
       CALL TRANSV(VRX,VTX,X1,Y1,UASYM,VASYM)
       UXY(I,J)= UASYM
       VXY(I,J)= VASYM
       ENDDO
       ENDDO
       ENDIF
C
C  Alternative--kept for testing the FT method
       DO J=1,NY
       DO I=1,NX
       X1= TDX(I)-XC
       Y1= TDY(J)-YC
       DIS= SQRT(X1*X1+Y1*Y1)
       IRR= INT(DIS/DR)
       IRR1= MAX0(IRR-2,1)
       IRR2= IRR1+MRF-1
       IRR2= MIN0(IRR2,MR)
       IRR2= MAX0(IRR2,MRF)
       IF (IRR2.GE.MR) IRR1=MR-MRF+1
       IF (IRR.LE.IRMN) THEN
       DO L=IRR1,IRR2
       XXR(L-IRR1+1)= FLOAT(L-1)*DR
       FFR(L-IRR1+1)= BVR(L)
       ENDDO 
       CALL SIINP(BVRX,DIS,FFR,XXR,MRF,0,0)
       DO L=IRR1,IRR2
       XXR(L-IRR1+1)= FLOAT(L-1)*DR
       FFR(L-IRR1+1)= BVT(L)
       ENDDO
       CALL SIINP(BVTX,DIS,FFR,XXR,MRF,0,0)
C--use BVRX,BVTX to get symmetric parts of U and V
       CALL TRANSV(BVRX,BVTX,X1,Y1,USYM,VSYM)
       if (IREMOVE.EQ.0) then
       UXY(I,J)= USYM
       VXY(I,J)= VSYM
       endif
       if (IREMOVE.EQ.1.AND.ISYM.NE.1) then
       UXY(I,J)= UXY(I,J) - USYM 
       VXY(I,J)= VXY(I,J) - VSYM
       endif
       ELSE
       UXY(I,J)= SPV
       VXY(I,J)= SPV
       ENDIF
       ENDDO
       ENDDO
C
C
       RETURN
       END

C-------------------------------------------------------------------
       SUBROUTINE ASYMS(QXY,NX,NY,TDX,TDY,XCS,YCS,IMMCH,SPV,IREMOVE)
       PARAMETER (MR=400,MT=361,MRF=6)  !'CHANGE ME'
       DOUBLE PRECISION XC,YC,DIS,XXR
       REAL SPV
       DIMENSION QXY(NX,NY),TDX(NX),TDY(NY)
       DIMENSION QRT(MR,MT),BUV(MR),BVR(MR),XXR(MRF),FFR(MRF)
C
C  finds the symmetric center of a scalar field and removes the
C  symmetric part
C
C  ON RETURN: 
C      QXY:  with the symmetric component removed and 
C            filled with SPV (except for zero) beyond MR
C  TDX, TDY: X,Y coordinates (monotonically increase)
C  XCS,YCS:  > 0 used as the coordinate without neibourhood searching; 
C            < 0 used as the grid index for neibourhood searching for
C            the center
C            -9999 for global searching of the center depending on IMMCH
C  IMMCH= 1 (0) searching maximum (minimum) defined as the center
C  Note:  Set MR=NX and MT=NY or so
C
       XC= XCS
       YC= YCS
       DR= ABS(TDX(2)-TDX(1))
C
c      print*,'ST'
       IF (SPV.EQ.0.) THEN
c      print*,'SPV = ',SPV
       DO J=1,NY
       DO I=1,NX
       SPV = SPV + QXY(I,J)
c      print*,'SPV = ',SPV
       ENDDO
       ENDDO
c      print*,'ST'
       SPV = SPV/(FLOAT(I)*FLOAT(J))
c       print *,'SPV= ',SPV
       ENDIF
C
       CALL CYLIN(QXY,NX,NY,TDX,TDY,QRT,MR,MT,BVR,DR,XC,YC,SPV,IMMCH)
       XCS=XC
       YCS=YC
c
c  the minimum radius for use
       TRX1= ABS(TDX(1)-XC)
       TRX2= ABS(TDX(NX)-XC)
       TRY1= ABS(TDY(1)-YC)
       TRY2= ABS(TDY(NY)-YC)
       TRMN= AMIN1(TRX1,TRX2,TRY1,TRY2)
       IRMN= INT(TRMN/DR)+1
c       print *,'Minimum radius dimension IRMN= ',IRMN
       IF (IRMN.GE.MR) THEN
c       print *,'!!  Need to increase MR...'
       ENDIF
C
C  Use the FT method to get wavenumber-N components
       ISYM = 1
       IF (ISYM.EQ.1.and.IREMOVE.NE.0) then
       DO I=1,MR
       CALL FT(QRT(I,1:MT-1),1,1,1,MT-1)
       ENDDO
       DO I=1,MR
       QRT(I,MT)= QRT(I,1)
       ENDDO
       CALL RTTOXY(QXY,NX,NY,TDX,TDY,QRT,MR,MT,BUV,DR,XC,YC,SPV,IMMCH)!
       ENDIF 
C
C  Alternative--kept for testing the FT method
       DO J=1,NY
         DO I=1,NX
           X1= TDX(I)-XC
           Y1= TDY(J)-YC
           DIS= SQRT(X1*X1+Y1*Y1)
           IRR= INT(DIS/DR)
           IRR1= MAX0(IRR-2,1)
           IRR2= IRR1+MRF-1
           IRR2= MIN0(IRR2,MR)
           IRR2= MAX0(IRR2,MRF)
           IF (IRR2.GE.MR) IRR1=MR-MRF+1
           IF (IRR.LE.IRMN) THEN
               DO L=IRR1,IRR2
                 XXR(L-IRR1+1)= FLOAT(L-1)*DR
                 FFR(L-IRR1+1)= BVR(L)
               ENDDO 
               DIS= DIS/DR
               DO II=1,MRF
                 XXR(II)= XXR(II)/DR
               ENDDO
               CALL SIINP(QR,DIS,FFR,XXR,MRF,0,0)
               if (IREMOVE.EQ.0) then
                 QXY(I,J)= QR
               endif
               if (IREMOVE.EQ.1.and.ISYM.NE.1) then
                 QXY(I,J)= QXY(I,J)- QR
               endif
           ELSE
             QXY(I,J)= SPV
           ENDIF
         ENDDO
       ENDDO
C
       SPV=0.
C
       RETURN
       END
C---------------------------------------------------------------------
       SUBROUTINE CYLIN(VXY,NX,NY,TDX,TDY,VRT,MR,MT,BVR,DR,
     &   XC,YC,SPV,IMMCH,ITER,TRT,HGT,HGT0,ITERXY,TER,
     &   ILS,XLAND,ALON,CLON)
       PARAMETER (INX=4,INY=4,NXF=5,NYF=5,TP=15)
       DOUBLE PRECISION XP,YP,XC,YC,XX,YY,XXF,YYF,PI,ANG
       DIMENSION VXY(NX,NY),VRT(MR,MT),TDX(NX),TDY(NY),BVR(MR)
       DIMENSION VV(INX,INY),XX(INX),YY(INY),TRT(MR,MT),TER(NX,NY)
       DIMENSION XXF(NXF),YYF(NYF),FF(NXF,NYF),TT(INX,INY)
       DIMENSION XLAND(NX,NY),ALON(NX,NY)
       REAL HGT,HGT0
       INTEGER ILS
       DATA PI/3.14159265358979/
C
C--PROGRAM TRANSFORM V(X,Y) TO V(R,THETA) for 2-D fields only.
C--INX,INY: total points for x and y interpolation.
C--NX,NY: VXY's dimensions of X and Y.
C--MR,MT: VRT's dimensions of R (radial) and THETA (tangential).
C  MT: for cyclic plot, MT is set to the partitions plus one 
C  MT= 5 means that there are five values around the cycle and
C  the first (J=1) and last values (J=MT) are the same.        
C--TDX,TDY: VXY's X and Y coordinates (rectangular).
C--XC,YC: the central point in TDX and TDY for the cylindrical origin.
C--DR: increment of radial distance (r=DR*(MR-1)).
C--ON RETURN: 
C     VRT(MR,MT): values at concentric points for analyses.
C     BVR(MR): the tangentially-averaged (symmetric) values
C
       NX1= NX-1
       NY1= NY-1
       MT1= MT-1
       IOUTSP= 0
       ICHECK= 0
       IMONO= 0
       NGRID= 20 
C       IF (NGRID.LE.0.01*NX.OR.NGRID.LE.0.01*NY) THEN
C         PRINT *,'!!  should consider larger values for NGRID= ', NGRID 
C       ENDIF
C
c--find xc,yc from max or min value
       IFIND= 0
       IF (XC.EQ.-9999.AND.YC.EQ.-9999.) THEN
         CALL MXMN2(VXY,NX,NY,1,NX,1,NY,IMX,JMY,IMMCH) 
         IFIND= 1
       ELSE
         IF (XC.LT.0.AND.YC.LT.0.) THEN
           INTXC= INT(ABS(XC))
           INTYC= INT(ABS(YC))
           IFIND= 1
         ENDIF
         IF (XC.GT.0.) THEN
           DO I=1,NX-1
C             print*,"I=,TDX(I)=,XC=,TDX(I+1)=",I,TDX(I),XC,TDX(I+1)
             IF (XC.GE.TDX(I).AND.XC.LT.TDX(I+1)) THEN
               INTXC= I
C              print*,"###############  I= ",I,"#################"
             ENDIF
           ENDDO
           IFIND= 0
         ENDIF
         IF (YC.GT.0.) THEN
           DO J=1,NY-1
             IF (YC.GE.TDY(J).AND.YC.LT.TDY(J+1)) THEN
               INTYC= J
             ENDIF
           ENDDO
           IFIND= 0
         ENDIF
         IF (INTXC.LE.0.OR.INTYC.LE.0) THEN
           PRINT *,'!!  Bad XC,YC...'
           RETURN
         ENDIF
         IA= MAX0(1,INTXC-NGRID)
         IB= MIN0(NX,INTXC+NGRID)
         JA= MAX0(1,INTYC-NGRID)
         JB= MIN0(NY,INTYC+NGRID)
         CALL MXMN2(VXY,NX,NY,IA,IB,JA,JB,IMX,JMY,IMMCH) 
       ENDIF
       PRINT *,'!!  interpolation center at imx,jmy: ',TDX(IMX),TDY(JMY)
       IF (IFIND.GE.1) THEN
         DO I=IMX-2,IMX+2
           XXF(I-IMX+3)= TDX(I)
         ENDDO
         DO J=JMY-2,JMY+2
           YYF(J-JMY+3)= TDY(J)
         ENDDO
         DO I=IMX-2,IMX+2
           DO J=JMY-2,JMY+2
             FF(I-IMX+3,J-JMY+3)= VXY(I,J)       
           ENDDO
         ENDDO
!!       print *,'XXF(i)=',XXF
!!       print *,'YYF(j)=',YYF
!!       DO J=1,NYF
!!       PRINT *,'FF(i,j)=',(FF(I,J),I=1,NXF)
!!       ENDDO
         CALL FINDMXY(NXF,NYF,XXF,YYF,FF,XC,YC,FC,IMMCH)
c         print *,'!!  finding new center coordinate XC,YC= ',XC,YC
       ENDIF
C
       DO II=1,NX
         IIP= MIN0(II+1,NX)
         IF (XC.GE.TDX(II).AND.XC.LT.TDX(IIP)) GO TO 88
       ENDDO
  88   IXC= II
       DO JJ=1,NY
         JJP= MIN0(JJ+1,NY)
         IF (YC.GE.TDY(JJ).AND.YC.LT.TDY(JJP)) GO TO 99
       ENDDO
  99   JYC= JJ
       IF (IXC.LE.2.OR.IXC.GE.NX1.OR.JYC.LE.2.OR.JYC.GE.NY1) THEN
         PRINT *,'!!  bad cylindrical origin at XC,YC!'
         RETURN
       ENDIF
C
C check XC,YC
       VXYR= VXY(IXC+1,JYC)-VXY(IXC,JYC)
       VXYL= VXY(IXC,JYC)-VXY(IXC-1,JYC)
       VXYU= VXY(IXC,JYC+1)-VXY(IXC,JYC)
       VXYD= VXY(IXC,JYC)-VXY(IXC,JYC-1)
       IF (VXYR*VXYL.GT.0.OR.VXYU*VXYD.GT.0.) THEN
         print *,'!!  XC,YC may not be the central position!'
       ENDIF
C
       print *,'!!  using final XC,YC= ',XC,YC
C
       DX= ABS(TDX(2)-TDX(1))
       DY= ABS(TDY(2)-TDY(1))
       DS= AMAX1(DX,DY)
       IF (DR.EQ.0.) DR=DS
       IF (DR.GE.10.*DS.OR.DR.LE.0.1*DS) THEN
         PRINT *,'!!  bad chosen radial increment!'
         RETURN
       ENDIF
C
C--first compute the first ring and others
C--representative of the central-point value for I=1
       DO 20 I=1,MR
         DO 10 J=1,MT
           RS= DR*FLOAT(I-1)
           ANG= 2.*PI/FLOAT(MT1)*FLOAT(J-1)  ! from zero degree
           XP= XC+RS*DCOS(ANG)  ! point on the ring
           YP= YC+RS*DSIN(ANG)  ! point on the ring
C--determine (I,J) of XP,YP on (X,Y)
           DO II=1,NX
             IIP= MIN0(II+1,NX)
             IF (XP.GE.TDX(II).AND.XP.LE.TDX(IIP)) GO TO 100
           ENDDO
           IXP= 9999
           GO TO 101
  100      IXP= II
  101      CONTINUE
           DO JJ=1,NY
             JJP= MIN0(JJ+1,NY)
             IF (YP.GE.TDY(JJ).AND.YP.LE.TDY(JJP)) GO TO 200
           ENDDO
           JYP= 9999
           GO TO 201
  200      JYP= JJ
  201      CONTINUE
C--extrapolation is avoided for exterior points
           IF (IXP.EQ.9999.OR.JYP.EQ.9999) THEN
             IF (IOUTSP.EQ.1) THEN
               PRINT *,'!!  Set special value at I,J= ',I,J
             ENDIF
             VRT(I,J)= SPV 
           ELSE 
             IXP1= IXP-INX/2+1
             IXP1= MAX0(IXP1,1)
             IXP2= IXP1+INX-1
             IF (IXP2.GT.NX) THEN
               IXP2= NX
               IXP1= IXP2-INX+1
             ENDIF
             DO II=IXP1,IXP2
               XX(II-IXP1+1)= TDX(II)
             ENDDO
             JYP1= JYP-INY/2+1
             JYP1= MAX0(JYP1,1)
             JYP2= JYP1+INY-1
             IF (JYP2.GT.NY) THEN
               JYP2= NY
               JYP1= JYP2-INY+1
             ENDIF
             DO JJ=JYP1,JYP2
               YY(JJ-JYP1+1)= TDY(JJ)
             ENDDO
             DO II=IXP1,IXP2
               DO JJ=JYP1,JYP2
                 VV(II-IXP1+1,JJ-JYP1+1)= VXY(II,JJ)
                 TT(II-IXP1+1,JJ-JYP1+1)= TER(II,JJ)
               ENDDO
             ENDDO
C--using bi-interpolation for interior points
c chi,202311: avoid the points under the terrain
             IF (ITERXY.EQ.1) THEN
               IF (ANY(HGT.LT.TT(:,:)).EQ..TRUE.) THEN 
                 VRT(I,J)= SPV
               ELSE
                 CALL BIINP(VP,XP,YP,VV,XX,YY,INX,INY,IMONO,0)
                 VRT(I,J)= VP
               ENDIF
               IF (ILS.EQ.1) THEN
                 IXP1= IXP-TP/2+1
                 IXP1= MAX0(IXP1,1)
                 IXP2= IXP1+TP-1
                 IF (IXP2.GT.NX) THEN
                   IXP2= NX
                   IXP1= IXP2-TP+1
                 ENDIF
                 JYP1= JYP-TP/2+1
                 JYP1= MAX0(JYP1,1)
                 JYP2= JYP1+TP-1
                 IF (JYP2.GT.NY) THEN
                   JYP2= NY
                   JYP1= JYP2-TP+1
                 ENDIF
                 IF (ANY(XLAND(IXP1:IXP2,JYP1:JYP2).LT.2.).EQ..TRUE.
     &               .and. ALON(IXP2,JYP2).LT.CLON) THEN
c                 print*,'!!!!! avoid the vicinity of Taiwan terrain'
                 VRT(I,J)= SPV
                 ENDIF
               ENDIF

             ELSE
               CALL BIINP(VP,XP,YP,VV,XX,YY,INX,INY,IMONO,0)
               VRT(I,J)= VP
             ENDIF
           ENDIF
 10      CONTINUE
 20    CONTINUE
c
c  the minimum radius dimension for use
c       TRX1= ABS(TDX(1)-XC)
c       TRX2= ABS(TDX(NX)-XC)
c       TRY1= ABS(TDY(1)-YC)
c       TRY2= ABS(TDY(NY)-YC)
c       TRMN= AMIN1(TRX1,TRX2,TRY1,TRY2)
c       IRMN= INT(TRMN/DR)+1
c       print *,'TRMN,IRMN= ',TRMN,IRMN
C
           IF (ICHECK.GE.2) THEN
             print *,' '
             print *,'The cylindrical values: '
             DO I=1,MR
               WRITE(6,*) (VRT(I,J),J=1,MT)
             ENDDO
           ENDIF
C
           DO I=2,MR
             TVR= 0.
             JSUM= 0.
             DO J=1,MT1
               IF (ITER.EQ.1) THEN
c                 IF (VRT(I,J).NE.SPV.AND.HGT.GE.TRT(I,J)) THEN
                 IM1= MAX0(I-1,1)
                 IP1= MIN0(I+1,MRSE)
                 JM1= MAX0(J-1,1)
                 JP1= MIN0(J+1,MT)
                 IF (HGT.GE.TRT(I,J).and.HGT.GE.TRT(IM1,J).and.
     &           HGT.GE.TRT(IP1,J).and.HGT.GE.TRT(I,JM1).and.
     &           HGT.GE.TRT(I,JP1).and.HGT0.GE.TRT(I,J).and.
     &           VRT(I,J).NE.SPV) THEN
                   JSUM= JSUM+1
                   TVR= TVR+ VRT(I,J)
                 ENDIF
               ELSE
                 IF (VRT(I,J).NE.SPV) THEN
                   JSUM= JSUM+1
                   TVR= TVR+ VRT(I,J)
                 ENDIF
               ENDIF
             ENDDO
             IF (JSUM.NE.0) THEN
               BVR(I)= TVR/FLOAT(JSUM)
             ELSE
               BVR(I)= SPV
             ENDIF
           ENDDO
           BVR(1)= VRT(1,1)
C
           IF (ICHECK.GE.1) THEN
             print *,' '
             print *,'The tangentially-averaged values: '
             WRITE(6,*) (BVR(I),I=1,MR)
           ENDIF
C
C
       RETURN
       END
C         
C---------------------------------------------------------------------
       SUBROUTINE RTTOXY(VXY,NX,NY,TDX,TDY,VRT,MR,MT,BVR,DR,
     &   XC,YC,SPV,IMMCH)
       PARAMETER (MRF=6,NXF=5,NYF=5)
       DOUBLE PRECISION XP,YP,XC,YC,XX,YY,XXF,YYF,PI,ANG,X1,Y1,DIS
       DIMENSION VXY(NX,NY),VRT(MR,MT),TDX(NX),TDY(NY),BVR(MR)
       DIMENSION VX1(MRF),VY1(MRF),XX(MRF),YY(MRF)
       DIMENSION XXF(NXF),YYF(NYF),FF(NXF,NYF)
       DATA PI/3.14159265358979/
C
C--PROGRAM TRANSFORM V(R,THETA) to V(X,Y) for 2-D fields only.
C--INX,INY: total points for x and y interpolation.
C--NX,NY: VXY's dimensions of X and Y.
C--MR,MT: VRT's dimensions of R (radial) and THETA (tangential).
C  MT: for cyclic plot, MT is set to the partitions plus one 
C  MT= 5 means that there are five values around the cycle and
C  the first (J=1) and last values (J=MT) are the same.        
C--TDX,TDY: VXY's X and Y coordinates (rectangular).
C--XC,YC: the central point in TDX and TDY for the cylindrical origin.
C--DR: increment of radial distance (r=DR*(MR-1)).
C--ON RETURN: 
C     VRT(MR,MT): values at concentric points for analyses.
C     BVR(MR): the tangentially-averaged (symmetric) values
C
       NX1= NX-1
       NY1= NY-1
       MT1= MT-1
       ICHECK= 0
       IMONO= 0
       NGRID= 20 
C
c--find xc,yc from max or min value
       IFIND= 0
       IF (XC.EQ.-9999.AND.YC.EQ.-9999.) THEN
       CALL MXMN2(VXY,NX,NY,1,NX,1,NY,IMX,JMY,IMMCH) 
       IFIND= 1
       ELSE
       IF (XC.LT.0.AND.YC.LT.0.) THEN
       INTXC= INT(ABS(XC))
       INTYC= INT(ABS(YC))
       IFIND= 1
       ENDIF 
       IF (XC.GT.0.) THEN
       DO I=1,NX-1
       IF (XC.GE.TDX(I).AND.XC.LT.TDX(I+1)) THEN
         INTXC= I
       ENDIF
       ENDDO
       IFIND= 0
       ENDIF
       IF (YC.GT.0.) THEN
       DO J=1,NY-1
       IF (YC.GE.TDY(J).AND.YC.LT.TDY(J+1)) THEN
         INTYC= J
       ENDIF
       ENDDO
       IFIND= 0
       ENDIF
       IA= MAX0(1,INTXC-NGRID)
       IB= MIN0(NX,INTXC+NGRID)
       JA= MAX0(1,INTYC-NGRID)
       JB= MIN0(NY,INTYC+NGRID)
       CALL MXMN2(VXY,NX,NY,IA,IB,JA,JB,IMX,JMY,IMMCH) 
       ENDIF
       IF (IFIND.GE.1) THEN
       DO I=IMX-2,IMX+2
       XXF(I-IMX+3)= TDX(I)
       ENDDO
       DO J=JMY-2,JMY+2
       YYF(J-JMY+3)= TDY(J)
       ENDDO
       DO I=IMX-2,IMX+2
       DO J=JMY-2,JMY+2
       FF(I-IMX+3,J-JMY+3)= VXY(I,J)       
       ENDDO
       ENDDO
       CALL FINDMXY(NXF,NYF,XXF,YYF,FF,XC,YC,FC,IMMCH)
       ENDIF
C
       DO II=1,NX
       IIP= MIN0(II+1,NX)
       IF (XC.GE.TDX(II).AND.XC.LT.TDX(IIP)) GO TO 88
       ENDDO
  88   IXC= II
       DO JJ=1,NY
       JJP= MIN0(JJ+1,NY)
       IF (YC.GE.TDY(JJ).AND.YC.LT.TDY(JJP)) GO TO 99
       ENDDO
  99   JYC= JJ
       IF (IXC.LE.2.OR.IXC.GE.NX1.OR.JYC.LE.2.OR.JYC.GE.NY1) THEN
       PRINT *,'!!  bad cylindrical origin at XC,YC!'
       RETURN
       ENDIF
C check XC,YC
       VXYR= VXY(IXC+1,JYC)-VXY(IXC,JYC)
       VXYL= VXY(IXC,JYC)-VXY(IXC-1,JYC)
       VXYU= VXY(IXC,JYC+1)-VXY(IXC,JYC)
       VXYD= VXY(IXC,JYC)-VXY(IXC,JYC-1)
       IF (VXYR*VXYL.GT.0.OR.VXYU*VXYD.GT.0.) THEN
       print *,'!!  XC,YC may not be the central position!'
       ENDIF
C
       print *,'In RTTOXY: XC,YC= ',XC,YC
C
       DX= ABS(TDX(2)-TDX(1))
       DY= ABS(TDY(2)-TDY(1))
       DS= AMAX1(DX,DY)
       IF (DR.EQ.0.) DR=DS
       IF (DR.GE.10.*DS.OR.DR.LE.0.1*DS) THEN
       PRINT *,'!!  bad chosen radial increment!'
       RETURN
       ENDIF
C
c  the minimum radius dimension for use
       TRX1= ABS(TDX(1)-XC)
       TRX2= ABS(TDX(NX)-XC)
       TRY1= ABS(TDY(1)-YC)
       TRY2= ABS(TDY(NY)-YC)
       TRMN= AMIN1(TRX1,TRX2,TRY1,TRY2)
       IRMN= INT(TRMN/DR)+1
       print *,'TRMN,IRMN= ',TRMN,IRMN
C
       IF (ICHECK.GE.2) THEN
       print *,' '
       print *,'The cylindrical values: '
       DO I=1,MR
       WRITE(6,*) (VRT(I,J),J=1,MT)
       ENDDO
       ENDIF
C
       DO I=2,MR
       TVR= 0.
       JSUM= 0.
       DO J=1,MT1
       IF (VRT(I,J).NE.SPV) THEN
       JSUM= JSUM+1
       TVR= TVR+ VRT(I,J)
       ENDIF
       ENDDO
       IF (JSUM.NE.0) THEN
         BVR(I)= TVR/FLOAT(JSUM)
       ELSE
         BVR(I)= SPV
       ENDIF
       ENDDO
       BVR(1)= VRT(1,1)
C
       IF (ICHECK.GE.1) THEN
       print *,' '
       print *,'The tangentially-averaged values: '
       WRITE(6,*) (BVR(I),I=1,MR)
       ENDIF
C
       DMT= 2.*PI/FLOAT(MT-1)
       DO 20 J=1,NY
       DO 10 I=1,NX
       X1= TDX(I)-XC
       Y1= TDY(J)-YC
       IF (ABS(Y1).LE.1.0E-3.AND.ABS(X1).LE.1.0E-3) THEN
       ANG= 0.
       VXY(I,J)= VRT(1,1)
       ELSE
       DIS= SQRT(X1*X1+Y1*Y1)
       ANG= DATAN2(Y1,X1)
       IF (ANG.LT.0.) ANG= 2.*PI + ANG 
       IRR= INT(DIS/DR)
       IRR1= MAX0(IRR-2,1)
       IRR2= IRR1+MRF-1
       IRR2= MIN0(IRR2,IRMN)
       IRR2= MAX0(IRR2,MRF)
       IF (IRR2.GE.IRMN) IRR1=IRR2-MRF+1
       IF (IRR1.EQ.1) IRR2=4
       IF (IRR.LE.IRMN) THEN
       JTT= INT(ANG/DMT)+1
       JTT1= MAX0(JTT-2,1)
       JTT2= JTT1+MRF-1
       DO II=IRR1,IRR2
       IF (JTT2.GE.MT) THEN
         DO JJ=JTT1,MT-1
         YY(JJ-JTT1+1)= FLOAT(JJ-1)*DMT
         VY1(JJ-JTT1+1)= VRT(II,JJ)  
         ENDDO
         JJ3= MT-JTT1+1
         DO JJ=JJ3,MRF
         YY(JJ)= YY(JJ3-1)+FLOAT(JJ-JJ3+1)*DMT
         VY1(JJ)= VRT(II,JJ-JJ3+1)  
         ENDDO
       ELSE 
         DO JJ=JTT1,JTT2
         YY(JJ-JTT1+1)= FLOAT(JJ-1)*DMT
         VY1(JJ-JTT1+1)= VRT(II,JJ)  
         ENDDO
       ENDIF
cc       print *,'JTT1,JTT2=',JTT1,JTT2
cc       print *,'ANG=',ANG
cc       print *,'DIS=',DIS
cc       print *,'VY1=',VY1
cc       print *,'YY=',YY
       CALL SIINP(VTX,ANG,VY1,YY,MRF,0,0)
       VX1(II-IRR1+1)= VTX
       XX(II-IRR1+1)= FLOAT(II-1)*DR
       ENDDO
       DIS= DIS/DR   ! better normalized for reducing round-off errors 
       DO II=1,MRF
       XX(II)= XX(II)/DR
       ENDDO
       IF (IRR1.LE.2) THEN 
         CALL SIINP(VRX,DIS,VX1,XX,4,0,0)
       ELSE
         CALL SIINP(VRX,DIS,VX1,XX,MRF,0,0)
       ENDIF
       VXY(I,J)= VRX
       ELSE
       VXY(I,J)= SPV
       ENDIF
       ENDIF
 10    CONTINUE
 20    CONTINUE
c
C
       RETURN
       END
C         
C------------------------------------------------------
       SUBROUTINE BIINP(QP,XP,YP,Q,X,Y,NX,NY,IMONO,IC)
       DOUBLE PRECISION XP,YP,X,Y,RX,RX0,RY,RY0
       DIMENSION Q(NX,NY),X(NX),Y(NY)
C--check coordinates which can be random without same points
C--kicking off is not performed for same (overlapped) points
       IF (IC.EQ.1) THEN
       DO II=1,NX
       DO JJ=1,NX
       IF (II.NE.JJ) THEN
       DX= X(II)-X(JJ)
       IF (ABS(DX).LE.1.0E-6) THEN
       PRINT *,'!!  Return due to ill coordinates in BIINP!'
       RETURN
       ENDIF
       ENDIF
       ENDDO
       ENDDO
       DO II=1,NY
       DO JJ=1,NY
       IF (II.NE.JJ) THEN
       DY= Y(II)-Y(JJ)
       IF (ABS(DY).LE.1.0E-6) THEN
       PRINT *,'!!  Return due to ill coordinates in BIINP!'
       RETURN
       ENDIF
       ENDIF
       ENDDO
       ENDDO
       ENDIF
C--using DO-LOOP rather than direct explicit expansion
       QP= 0.
       DO 20 II=1,NX
       DO 10 JJ=1,NY
       RX0= 1.
       DO I=1,NX
       IF (II.NE.I) THEN
       RX= (XP-X(I))/(X(II)-X(I))*RX0
       RX0= RX
       ENDIF
       ENDDO
       RY0= 1.
       DO J=1,NY
       IF (JJ.NE.J) THEN
       RY= (YP-Y(J))/(Y(JJ)-Y(J))*RY0
       RY0= RY
       ENDIF
       ENDDO
       QP= QP + RX*RY*Q(II,JJ)
 10    CONTINUE
 20    CONTINUE
C
C
       if (imono.eq.1) then
       QMIN= Q(1,1)
       QMAX= Q(1,1)
       DO II=1,NX
       DO JJ=1,NY
       QMIN= AMIN1(QMIN,Q(II,JJ))
       QMAX= AMAX1(QMAX,Q(II,JJ))
       QP= AMAX1(QP,QMIN)
       QP= AMIN1(QP,QMAX)
       ENDDO
       ENDDO
       endif
C
C
       RETURN
       END
C
C------------------------------------------------------
       SUBROUTINE SIINP(QP,XP,Q,X,NX,IMONO,IC)
       DOUBLE PRECISION XP,X,RX,RX0
       DIMENSION Q(NX),X(NX)
C--check coordinates which can be random without same points
C--kicking off is not performed for same (overlapped) points
       IF (IC.EQ.1) THEN
       DO II=1,NX
       DO JJ=1,NX
       IF (II.NE.JJ) THEN
       DX= X(II)-X(JJ)
       IF (ABS(DX).LE.1.0E-6) THEN
       PRINT *,'!!  Return due to ill coordinates in BIINP!'
       RETURN
       ENDIF
       ENDIF
       ENDDO
       ENDDO
       ENDIF
C--using DO-LOOP rather than direct explicit expansion
       QP= 0.
       DO 20 II=1,NX
       RX0= 1.
       DO I=1,NX
       IF (II.NE.I) THEN
       RX= (XP-X(I))/(X(II)-X(I))*RX0
       RX0= RX
       ENDIF
       ENDDO
       QP= QP + RX*Q(II)
 20    CONTINUE
C
C
       if (imono.eq.1) then
       QMIN= Q(1)
       QMAX= Q(1)
       DO II=1,NX
       QMIN= AMIN1(QMIN,Q(II))
       QMAX= AMAX1(QMAX,Q(II))
       QP= AMAX1(QP,QMIN)
       QP= AMIN1(QP,QMAX)
       ENDDO
       endif
C
C
       RETURN
       END
C
C--------------------------------------------------------------
        SUBROUTINE BILIN(VP,XP,YP,V11,V21,V12,V22,X1,X2,Y1,Y2,IC)       
C--simple direct expansion used for boundary points if necessary
       IF (IC.EQ.1) THEN
       IF (X1.EQ.X2.OR.Y1.EQ.Y2) THEN
       PRINT *,'!!  Return due to ill coordinates in BILIN!'
       RETURN
       ENDIF
       ENDIF
       X12= X1-X2
       Y12= Y1-Y2
       VP= (XP-X2)*(YP-Y2)/(X12*Y12)*V11+
     &     (XP-X1)*(YP-Y2)/(-X12*Y12)*V21+
     &     (XP-X2)*(YP-Y1)/(-X12*Y12)*V12+
     &     (XP-X1)*(YP-Y1)/(X12*Y12)*V22   
C
C
       RETURN
       END
C
C------------------------------------------------------------------
       FUNCTION CUBLA3(XP,NX,TDX,Q,CHECK)
       DIMENSION TDX(NX),Q(NX)
C
C  for stretched grids with unknown XP
c
       CHECK= 0.
C
       N= NX
       X= XP
C
C  determine I for X
C
       IF (X.LE.TDX(1)) THEN
         IIP= 1
         GO TO 10
       ENDIF
       IF (X.GE.TDX(N)) THEN
         IIP= N-1
         GO TO 10
       ENDIF
       DO II=1,N-1
       IF (X.GE.TDX(II).AND.X.LT.TDX(II+1)) THEN
         IIP= II
         GO TO 10
       ENDIF
       ENDDO
C
  10   IIP= MIN0(IIP,N)
C
       I= IIP
C
       IF (I.GE.2.AND.I.LE.N-2) THEN
C
       IM1= MAX0(I-1,1)
       IP1= MIN0(I+1,N)
       IP2= MIN0(I+2,N)
       X0= TDX(IM1)
       X1= TDX(I)
       X2= TDX(IP1)
       X3= TDX(IP2)
C
       XX1= X - X1
       XX2= X - X2
       XX3= X - X3
       XX0= X - X0
       X10= X1 - X0
       X20= X2 - X0
       X30= X3 - X0
       X21= X2 - X1
       X32= X3 - X2
       X31= X3 - X1
C
       P3X= -(XX1/X10)*(XX2/X20)*(XX3/X30)*Q(IM1)
     $   +  (XX0/X10)*(XX2/X21)*(XX3/X31)*Q(I)
     $   - (XX0/X20)*(XX1/X21)*(XX3/X32)*Q(IP1)
     $   + (XX0/X30)*(XX1/X31)*(XX2/X32)*Q(IP2)
C
       ELSE

       IP1= MAX0(I+1,1)
       P3X= ( (X-TDX(I))*Q(IP1) + (TDX(IP1)-X)*Q(I) )/(TDX(IP1)-TDX(I))

       ENDIF

       CUBLA3= P3X
C
       RETURN
       END
C
C---------------------------------------------------------------------
       SUBROUTINE SVRVT(U,V,VR,VT,MR,MT,MZ)
       DOUBLE PRECISION PI,ANG,ANG1,ANG2,ANG3,ANG4
       DIMENSION U(MR,MT,MZ),V(MR,MT,MZ),VR(MR,MT,MZ),VT(MR,MT,MZ)
       DATA PI/3.14159265358979/
C
C--this program decompose the 2-D wind to radial and tangential
Ccomponents
C  U,V: the 2-D wind on the cylindrical coordinate.
C  MR,MT: the dimensions in radial and tangential directions.
C  MZ:    the vertical dimension (set to 1 for 2-D fields).
C  ON RETUEN:
C     VR,VT: the radial and tangential components.
C
C     VR=  U*cos(a)+V*sin(a)
C     VT= -U*sin(a)+V*cos(a)
C
C
       icheck= 0
       ischeme= 2
C
       IF (ISCHEME.EQ.1) THEN     !  illustrative form
       DO K=1,MZ
       DO I=1,MR
       DO J=1,MT
c chi,202311: avoid the -999 points
       IF (U(I,J,K).EQ.-999.) THEN
         VR(I,J,K)= -999.
         VT(I,J,K)= -999.
c chi,202311: end
       ELSE 
       ANG= 2.*PI/FLOAT(MT-1)*FLOAT(J-1)  ! from zero degree
       IF (ANG.GE.0.AND.ANG.LE.0.5*PI) THEN
       ANG1= ANG
       UUR= U(I,J,K)*DCOS(ANG1)
       UUT= -U(I,J,K)*DSIN(ANG1)
       VVR= V(I,J,K)*DSIN(ANG1)
       VVT= V(I,J,K)*DCOS(ANG1)
       ENDIF
       IF (ANG.GT.0.5*PI.AND.ANG.LE.PI) THEN
       ANG2= PI-ANG
       UUR= -U(I,J,K)*DCOS(ANG2)
       UUT= -U(I,J,K)*DSIN(ANG2)
       VVR= V(I,J,K)*DSIN(ANG2)
       VVT= -V(I,J,K)*DCOS(ANG2)
       ENDIF
       IF (ANG.GT.PI.AND.ANG.LE.1.5*PI) THEN
       ANG3= 1.5*PI-ANG
       UUR= -U(I,J,K)*DSIN(ANG3)
       UUT= U(I,J,K)*DCOS(ANG3)
       VVR= -V(I,J,K)*DCOS(ANG3)
       VVT= -V(I,J,k)*DSIN(ANG3)
       ENDIF
       IF (ANG.GT.1.5*PI.AND.ANG.LE.2.*PI) THEN
       ANG4= 2.*PI-ANG
       UUR= U(I,J,K)*DCOS(ANG4)
       UUT= U(I,J,K)*DSIN(ANG4)
       VVR= -V(I,J,K)*DSIN(ANG4)
       VVT= V(I,J,K)*DCOS(ANG4)
       ENDIF
       VR(I,J,K)= UUR+VVR
       VT(I,J,K)= UUT+VVT
       ENDIF
       ENDDO
       ENDDO
       ENDDO
       ENDIF
C
       IF (ISCHEME.EQ.2) THEN      ! compact form
       DO K=1,MZ
       DO I=1,MR
       DO J=1,MT
c chi,202311: avoid the -999 points
       IF (U(I,J,K).EQ.-999.) THEN
         VR(I,J,K)= -999.
         VT(I,J,K)= -999.
c chi,202311: end
       ELSE
       ANG= 2.*PI/FLOAT(MT-1)*FLOAT(J-1)  ! from zero degree
       VR(I,J,K)=  U(I,J,K)*dcos(ang)+V(I,J,K)*dsin(ang)
       VT(I,J,K)= -U(I,J,K)*dsin(ang)+V(I,J,K)*dcos(ang)
       ENDIF
       ENDDO
       ENDDO
       ENDDO
       ENDIF
C
       if (icheck.eq.1) then
       print *,'VR='
       DO I=1,MR
       WRITE(6,*) (VR(I,J,1),J=1,MT)
       ENDDO
       print *,'VT='
       DO I=1,MR
       WRITE(6,*) (VT(I,J,1),J=1,MT)
       ENDDO
       endif
C
C
       RETURN
       END
C---------------------------------------------------------------------
       subroutine findmxy(nx,ny,x,y,f,xp,yp,fp,immch)
c  program to find the position of extreme values with given points
       parameter (np=10)
       double precision xp,yp,x,y
       dimension x(nx),y(ny),f(nx,ny)
c
       if (immch.eq.0) then
         ifmax= 0
         ifmin= 1
       else
         ifmax= 1
         ifmin= 0
       endif
       imono= 0
       icheck= 1
       errf= 1.0e-10
c
       print *,'!!  zooming with nx,ny= ',nx,ny,' with step np= ',np
c
       if (icheck.eq.1) then
       icmin= 0
       icmax= 0
       fmin= f(1,1)  
       fmax= f(1,1)  
       do i=1,nx
       do j=1,ny
       fmin= amin1(f(i,j),fmin)
       fmax= amax1(f(i,j),fmax)
       enddo
       enddo
       do i=1,nx
       do j=1,ny
       if (f(i,j).le.fmin) then
         imin= i
         jmin= j
         icmin= icmin + 1
       endif
       if (f(i,j).ge.fmax) then
         imax= i
         jmax= j
         icmax= icmax + 1
       endif
       enddo
       enddo
       if (icmin.gt.1) then
       print *,'!!  min. value has no single point...'
       else
       print *,'!!  min. grid value= ',fmin,' at ii,jj=',imin,jmin
       endif
       if (icmax.gt.1) then
       print *,'!!  max. value has no single point...'
       else
       print *,'!!  max. grid value= ',fmax,' at ii,jj=',imax,jmax
       endif
       endif
c
c  find position of the real extreme value using searching stepping
       if (ifmin.eq.1) then 
       icc= 0
       xpmin= x(imin) 
       ypmin= y(jmin)
       fmin= f(imin,jmin)
       imin1= max0(imin-1,1)
       imin2= min0(imin+1,nx)
       jmin1= max0(jmin-1,1)
       jmin2= min0(jmin+1,ny)
       x2= x(imin2)
       x1= x(imin1)
       y2= y(jmin2)
       y1= y(jmin1)
 10    fmino= fmin
       xss= (x2-x1)/float(np)
       yss= (y2-y1)/float(np)
       do n=1,np
       do m=1,np
       xp= x1+ xss*float(n)
       yp= y1+ yss*float(m)
       call biinp(fp,xp,yp,f,x,y,nx,ny,imono,ic)
!!       print *,'fp,fmin=',fp,fmin
       if (fp.le.fmin) then
         xpmin= xp
         ypmin= yp
         fmin= fp
       endif 
       enddo
       enddo
       if (abs(fmin-fmino).ge.errf) then
         icc= icc +1
         print *,'!!  fmin,fmino= ',fmin,fmino,' at times=',icc
         x1= xpmin-xss
         x2= xpmin+xss
         y1= ypmin-yss
         y2= ypmin+yss
         go to 10
       endif
       print *,'!!  found min. value= ',fmin,' at x,y=',xpmin,ypmin
       xp= xpmin
       yp= ypmin
       endif
c
       if (ifmax.eq.1) then 
       jcc= 0
       xpmax= x(imax)
       ypmax= y(jmax)
       fmax= f(imax,jmax)
       imax1= max0(imax-1,1)
       imax2= min0(imax+1,nx)
       jmax1= max0(jmax-1,1)
       jmax2= min0(jmax+1,ny)
       x2= x(imax2)
       x1= x(imax1)
       y2= y(jmax2)
       y1= y(jmax1)
 20    fmaxo= fmax
       xss= (x2-x1)/float(np)
       yss= (y2-y1)/float(np)
       do n=1,np
       do m=1,np
       xp= x1+ xss*float(n)
       yp= y1+ yss*float(m)
       call biinp(fp,xp,yp,f,x,y,nx,ny,imono,ic)
       if (fp.ge.fmax) then
         xpmax= xp
         ypmax= yp
         fmax= fp
       endif 
       if (abs(fmax-fmaxo).ge.errf) then
         jcc= jcc+ 1
         print *,'!!  fmax,fmaxo= ',fmax,fmaxo,' at times=',jcc
         x1= xpmax-xss
         x2= xpmax+xss
         y1= ypmax-yss
         y2= ypmax+yss
         go to 20
       endif
       enddo
       enddo
       print *,'!!  found max. value= ',fmax,' at x,y=',xpmax,ypmax
       xp= xpmax
       yp= ypmax
       endif
c
c
       return 
       end
C
C--------------------------------------------------------------------
       SUBROUTINE MXMN2(F,NX,NY,IA,IB,JA,JB,IMX,JMY,IMMCH)
C  routine find max and min on grids
       DIMENSION F(NX,NY)
C
       IF (IA.LT.1.OR.IB.GT.NX.OR.JA.LT.1.OR.JB.GT.NY) THEN
       print *,'!!  improper domain range chosen...'
       stop
       ENDIF
C
       icmin= 0
       icmax= 0
       fmin= f(IA,JA)
       fmax= f(IA,JA)
       do i=ia,ib
       do j=ja,jb
       fmin= amin1(f(i,j),fmin)
       fmax= amax1(f(i,j),fmax)
       enddo
       enddo
       do i=ia,ib
       do j=ja,jb
       if (f(i,j).eq.fmin) then
         imin= i
         jmin= j
         icmin= icmin + 1
       endif
       if (f(i,j).eq.fmax) then
         imax= i
         jmax= j
         icmax= icmax + 1
       endif
       enddo
       enddo
       print *,'!!  searching in the zone of IA,IB,JA,JB: ',IA,IB,JA,JB
       if (icmin.gt.1) then
       print *,'!!  min. value has no single point...'
       else
       print *,'!!  min. grid value= ',fmin,' at i,j=',imin,jmin
       endif
       if (icmax.gt.1) then
       print *,'!!  max. value has no single point...'
       else
       print *,'!!  max. grid value= ',fmax,' at i,j=',imax,jmax
       endif
c
       if (immch.eq.1) then
         imx= imax
         jmy= jmax
       endif
c
       if (immch.eq.0) then
         imx= imin
         jmy= jmin
       endif
C
       print *,'!!  returning imx,jmy= ',imx,jmy
c
       RETURN
       END

C--------------------------------------------------------------------
       SUBROUTINE MXMN3(F,NX,NY,NZ,IMX,JMY,KMZ,IMMCH)
C  routine find max and min on grids
       DIMENSION F(NX,NY,NZ)
C
       icmin= 0
       icmax= 0
       fmin= f(1,1,1)
       fmax= f(1,1,1)
       do i=1,nx
       do j=1,ny
       do k=1,nz
       fmin= amin1(f(i,j,k),fmin)
       fmax= amax1(f(i,j,k),fmax)
       enddo
       enddo
       enddo
       do i=1,nx
       do j=1,ny
       do k=1,nz
       if (f(i,j,k).eq.fmin) then
         imin= i
         jmin= j
         kmin= k
         icmin= icmin + 1
       endif
       if (f(i,j,k).eq.fmax) then
         imax= i
         jmax= j
         kmax= k
         icmax= icmax + 1
       endif
       enddo
       enddo
       enddo
       if (icmin.gt.1) then
       print *,'!!  min. value has no single point...'
       else
       print *,'!!  min. grid value= ',fmin,' at i,j,k=',imin,jmin,kmin
       endif
       if (icmax.gt.1) then
       print *,'!!  max. value has no single point...'
       else
       print *,'!!  max. grid value= ',fmax,' at i,j,k=',imax,jmax,kmax
       endif
c
       if (immch.eq.1) then
         imx= imax
         jmy= jmax
         kmz= kmax
       endif
c
       if (immch.eq.0) then
         imx= imin
         jmy= jmin
         kmz= kmin
       endif
C
       print *,'!!  returning imx,jmy,kmz= ',imx,jmy,kmz
c
       RETURN
       END

C----------------------------------------------------------
       SUBROUTINE TRANSV(VR,VT,X,Y,U,V)
       DOUBLE PRECISION X,Y,ANG
C--transfer VR,VT back to U,V with a given location 
C  U= -VT*sin(a) + VR*cos(a)
C  V=  VT*cos(a) + VR*sin(a)
C
       IF (Y.EQ.0.AND.X.EQ.0.) THEN
       ANG= 0.
       ELSE
       ANG= DATAN2(Y,X)
       ENDIF
C
       U= -VT*SIN(ANG)+VR*COS(ANG)
       V=  VT*COS(ANG)+VR*SIN(ANG)
C
C
       RETURN
       END
C
c--------------------------------------------------------------------
c Discrete Fourier Transform
c Modifying from the code of fft by "Rosetta code"
      subroutine FT(xx,w,m,d,MT)
      integer,       parameter :: dp=selected_real_kind(15,300)
      real,parameter    :: pi=3.141592653589793238460_dp
      complex(kind=dp),dimension(MT)  :: x
      complex(kind=dp),dimension(MT)  :: x1
      integer                                        :: MT,i,w,m,d,ps
      real :: xx(MT)
c
c FT(xx,w,m,d,MT)
c program for discrete Fourier transform and
c inverse discrete Fourier transform
c xx   : input/output matrix
c w   : wno. low filter(wno. filter)
c m   : wno. high filter
c d   : working type
c MT   : the length of your matrix
c
c d=1 : only reserving the signal of wave number "w" in Inverse DFT
c d=2 : reserving the signal from wave number "w" to "m"
c d=3 : only reserving the signal of wave "w" and "m"
c
      ps=0
      x=cmplx(xx)
      call fft(x,MT)
      do i=1,MT
       x1(i)=0.
      enddo
c      
      if(d .eq. 1)then
       if(w .eq. 0)then
        x1(w+1)=x(w+1)
       else 
        x1(w+1)=x(w+1)
        x1(MT-w+1)=x(MT-w+1)
       endif
      elseif(d .eq. 2)then
       if(w .eq. 0)then
        x1(2+w-1:2+m-1)=x(2+w-1:2+m-1)
        x1(MT-m+1:MT)=x(MT-m+1:MT)
       else
        x1(2+w-1:2+m-1)=x(2+w-1:2+m-1)
        x1(MT-m+1:MT-w+1)=x(MT-m+1:MT-w+1)
       endif
      elseif(d .eq. 3)then
       if(w .eq. 0)then
        x1(w+1)=x(w+1)
        x1(m+1)=x(m+1)
        x1(MT-m+1)=x(MT-m+1)
       else
        x1(w+1)=x(w+1)
        x1(MT-w+1)=x(MT-w+1)
        x1(m+1)=x(m+1)
        x1(MT-m+1)=x(MT-m+1)
       endif
      elseif(d .eq. 0)then
       print*,'test...ps'
       x1(:)=x(:)
       go to 82
      else
       print*,'something  wrong on fft!'
       print*,'the variable d should be 0,1 or 2 or 3'
       return
      endif
  82  CONTINUE
       if(ps .eq. 1)then
c  82  CONTINUE
       do i=1,MT
        x1(i)=x1(i)*conjg(x1(i))/MT
       enddo
       go to 74
       endif
      call ifft(x1,MT)
  74  CONTINUE
      do i=1,MT
       xx(i)=x1(i)/MT
      enddo
      return
      end subroutine FT
C
c----------------------------------------------------------------------
      subroutine fft(x,MT)
      implicit none
      integer,       parameter :: dp=selected_real_kind(15,300)
      real,parameter    :: pi=3.141592653589793238460_dp
      complex(kind=dp), dimension(MT), intent(inout)  :: x
      complex(kind=dp)                               :: t
      integer                                        :: MT
      integer                                        :: i,j
      complex(kind=dp), dimension(MT)    :: y
c
      do i=1,MT
       y(i)=0.
       do j=1,MT
        t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(j-1,dp)*real(i-1,dp)
     $    /real(MT,dp),kind=dp))*x(j)
        y(i)     = y(i) + t
       enddo
      enddo
c
      do i=1,MT
       x(i)=y(i)
      enddo
c
      return
      end subroutine fft
c
c--------------------------------------------------------------------
      subroutine ifft(x,MT)
      implicit none
      integer,       parameter :: dp=selected_real_kind(15,300)
      real,parameter    :: pi=3.141592653589793238460_dp
      complex(kind=dp), dimension(MT), intent(inout)  :: x
      complex(kind=dp)                               :: t
      integer                                        :: MT
      integer                                        :: i,j
      complex(kind=dp), dimension(MT)    :: y

      if(MT .le. 1) return

!
      do i=1,MT
       y(i)=0.
       do j=1,MT
        t=exp(cmplx(0.0_dp,2.0_dp*pi*real(j-1,dp)*real(i-1,dp)
     $    /real(MT,dp),kind=dp))*x(j)
        y(i)     = y(i) + t
       enddo
      enddo
!
      do i=1,MT
       x(i)=real(y(i))
      enddo

      return
      end subroutine ifft

C
C--------------------------------------------------------------------
       SUBROUTINE MAXMIN(FIELD,VAR,MX,MY,MZ1,MZ2,NX,NY,NZ,IBC,ICHK)
       CHARACTER FIELD*(*)
       DIMENSION VAR(NX,NY,NZ)
C
       VMIN= 1.0E+29
       VMAX= -1.0E+29
C
       IF (IBC.EQ.1) then
       I1= 2
       I2= MX-1
       J1= 2
       J2= MY-1
       ELSE
       I1= 1
       I2= MX
       J1= 1
       J2= MY
       ENDIF
       
       DO K=MZ1,MZ2
       DO J=J1,J2
       DO I=I1,I2
       VMIN= AMIN1(VMIN,VAR(I,J,K))
       VMAX= AMAX1(VMAX,VAR(I,J,K))
       ENDDO
       ENDDO
       ENDDO
C
       IF (ICHK.EQ.1) THEN
       DO K=MZ1,MZ2
       DO J=J1,J2
       DO I=I1,I2
       IF (VAR(I,J,K).EQ.VMIN) THEN
         IMIN= I
         JMIN= J
         KMIN= K
       ENDIF
       IF (VAR(I,J,K).EQ.VMAX) THEN
         IMAX= I
         JMAX= J
         KMAX= K
       ENDIF
       ENDDO
       ENDDO
       ENDDO
       ENDIF
C
       PRINT *,FIELD,' from ',VMIN,' to ',VMAX
C
       IF (ICHK.EQ.1) THEN
       WRITE(6,1001) IMIN,JMIN,KMIN,IMAX,JMAX,KMAX,MZ1,MZ2
 1001  FORMAT(5X,'>>> Min at ',I3,1X,I3,1X,I3,
     &        2X,'&  Max at ',I3,1X,I3,1X,I3,'  for K=',I3,' to ',I3)
       ENDIF
C
C
       RETURN
       END
C
CC--------------------------------------------------------------------
       SUBROUTINE MAXMIN2(FIELD,VAR,MX,MZ1,MZ2,NX,NZ,IBC,ICHK)
       CHARACTER FIELD*(*)
       DIMENSION VAR(NX,NZ)
C
       IF (IBC.EQ.1) then
       I1= 2
       I2= MX-1
       ELSE
       I1= 1
       I2= MX
       ENDIF
C
       VMIN= 1.0E+29
       VMAX= -1.0E+29
       DO K=MZ1,MZ2
       DO I=I1,I2
       VMIN= AMIN1(VMIN,VAR(I,K))
       VMAX= AMAX1(VMAX,VAR(I,K))
       ENDDO
       ENDDO
C
       IF (ICHK.EQ.1) THEN
       DO K=MZ1,MZ2
       DO I=I1,I2
       IF (VAR(I,K).EQ.VMIN) THEN
         IMIN= I
         KMIN= K
       ENDIF
       IF (VAR(I,K).EQ.VMAX) THEN
         IMAX= I
         KMAX= K
       ENDIF
       ENDDO
       ENDDO
       ENDIF
C
       PRINT *,FIELD,' from ',VMIN,' to ',VMAX
C
       IF (ICHK.EQ.1) THEN
       WRITE(6,1001) IMIN,KMIN,IMAX,KMAX,MZ1,MZ2
 1001  FORMAT(5X,'>>> Min at ',I3,1X,I3,
     &        2X,'&  Max at ',I3,1X,I3,'  for K=',I3,' to ',I3)
       ENDIF
C
C
       RETURN
       END
C
