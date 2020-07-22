#include "TIDES_OPTIONS.h"

C TODO:
C 1. Replace call to sunLon and sunLat to avoid using SPICE
C 2. Replace GM using Kepler's third law 
C 3. Add eccentricity and obliquity as input parameters
C 4. Modify calculation for eccentricity and obliquity tide
C 5. Decide on how to add the k-loop? Probably multiply forcing by
C    scaling factor. 

CBOP
C !ROUTINE: TIDES_CALC

C !INTERFACE: ==========================================================
      SUBROUTINE TIDES_CALC(
     I                     myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate custom tendency terms outside k-loop in DO_OCEANIC_PHYS

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "TIDES_PARAMS.h"

C !INPUT PARAMETERS: ===================================================
C myTime    :: Current time in simulation
C myIter    :: Current time-step number
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================
CEOP

#ifdef ALLOW_TIDES

C !FUNCTIONS: ==========================================================
      EXTERNAL MUZEN
      REAL*8   MUZEN
      EXTERNAL LEGENDRE2
      _RL      LEGENDRE2

C !LOCAL VARIABLES: ====================================================
C bi,bj     :: Tile indices
      INTEGER bi, bj, i, j
      INTEGER myDate(4)
      INTEGER yr, mon, day, ss, ld, wd, H, M, S
      CHARACTER*20 timestr
      REAL*8 time, posrec(3), rkm, lt, flat, alt, rng, rasc, dec
      REAL*8 lon, lat
      REAL*8 sinLat, cosLat, cosM, sinM, cos2Lon, sin2Lon, sin2Lat
      REAL*8 cosLon
      REAL*8 eccMag, oblMag
      _RL muSU, muMO, rSU, rMO, xi, VSU, VMO, Vperm, cent

      _BEGIN_MASTER(myThid)

      flat = 0D0
      rkm = tides_earthRadius/1 _d 3

C     CALL cal_GetDate( myIter, myTime, myDate, myThid )
      yr = myDate(1)/10000
      mon = MOD(myDate(1)/100,100)
      day = MOD(myDate(1),100)
      H = myDate(2)/10000
      M = MOD(myDate(2)/100,100)
      S = MOD(myDate(2),100)
C     convert UTC to ephemeris time using SPICE
C      WRITE(timestr,'(I4,X,I2.2,X,I2.2,X,I2.2,X,I2.2,X,I2.2)'),
C     &   yr, mon, day, H, M, S
C     CALL STR2ET(timestr, time)

C     CALL SPKPOS('SUN',time,tides_orientData,'NONE','EARTH',posrec,lt)
C     CALL RECGEO(posrec,rkm,flat,lonSU,latSU,alt)
C     CALL RECRAD(posrec,rng,rasc,dec)
C     convert to m
      rSU = rng*1D3
      xi = tides_earthRadius/rSU
      aSU = -tides_sunGM/rSU*(xi*xi)

C     CALL SPKPOS('MOON',time,tides_orientData,'NONE','EARTH',posrec,lt)
C     CALL RECGEO(posrec,rkm,flat,lonMO,latMO,alt)
C     CALL RECRAD(posrec,rng,rasc,dec)
C     convert to m
      rMO = rng*1D3
      xi = tides_earthRadius/rMO
      aMO = -tides_moonGM/rMO*(xi*xi)

C     time is already relative to J2000 epoch, convert to centuries
      cent = time/36525 _d 0/86400 _d 0
      tides_permC = tides_permC0 + tides_permC1*cent

C      _END_MASTER(myThid)
C      _BARRIER

C      DO bj=myByLo(myThid),myByHi(myThid)
C       DO bi=myBxLo(myThid),myBxHi(myThid)
C        DO j=0,sNy+1
C         DO i=0,sNx+1
C          lon = XC(i,j,bi,bj)*deg2rad
C          lat = YC(i,j,bi,bj)*deg2rad
C          muSU = MUZEN(lon, lat, lonSU, latSU)
C          muMO = MUZEN(lon, lat, lonMO, latMO)
C          VSU = aSU*LEGENDRE2(muSU)
C          VMO = aMO*LEGENDRE2(muMO)
C          Vperm = tides_permC*SQRT(5 _d 0)*LEGENDRE2(SIN(lat))
C          phiTide2d(i,j,bi,bj) = (1 _d 0 + tides_k2 - tides_h2)
C     &                          *(VSU + VMO - Vperm)
C         ENDDO
C        ENDDO
C       ENDDO
C      ENDDO

C     Calculate tidal potential (to first order) in 
C     eccentricity and obliquity.
      cosM = 0 _d 0
      sinM = 0 _d 0
   
      eccMag = 0 _d 75 * tides_earthRadius*tides_earthRadius
     &         *omega*omega*tides_ecc
      oblMag = -1 _d 5 * tides_earthRadius*tides_earthRadius
     &          *omega*omega*tides_obl*deg2rad

      _END_MASTER(myThid)
      _BARRIER


      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=0,sNy+1
         DO i=0,sNx+1
          lon = XC(i,j,bi,bj)*deg2rad
          lat = YC(i,j,bi,bj)*deg2rad
          
          sinLat  = SIN(lat)
          cosLat  = COS(lat)
          cosLon  = COS(lon)
          cos2Lon = COS(2 _d 0 * lon)
          sin2Lon = SIN(2 _d 0 * lon)
          sin2Lat = SIN(2 _d 0 * lat)

          phiTide2d(i,j,bi,bj) = (1 _d 0 + tides_k2 - tides_h2)
C               eccentricity forcing
     &               *( ( (1 _d 0 - 3 _d 0 * sinLat*sinLat)*cosM
     &                    + cosLat*cosLat * (3 _d 0 *cosM*cos2Lon 
     &                                     + 4 _d 0 *sinM*sin2Lon)
     &                   )*eccMag
C               obliquity forcing
     &                 + (cosM * sin2Lat * cosLon)*oblMag )
      
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_TIDES */

      RETURN
      END

      FUNCTION MUZEN(lon1, lat1, lon2, lat2)
      REAL*8 lon1, lat1, lon2, lat2, muzen
      MUZEN = SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1)
      RETURN
      END

      FUNCTION LEGENDRE2(mu)
      _RL LEGENDRE2, mu
      LEGENDRE2 = 1.5 _d 0*mu*mu - 0.5 _d 0
      RETURN
      END
