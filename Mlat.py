from GPS_utils import GPS_utils
import math

#4 	Frankfurt 	 	yes 	50.04854 	8.48795 	96.1 
#18 	Zurich 	 	yes 	47.45382 	8.55791 	470
#19 	Basel-Mulhouse 		yes 	47.59975 	7.53157 	314.7 

#"1623924000004404293:4,1623924000003618910:18,1623924000003779350:19",3,47.684189,8.524155
Receivers=[[50.04854,8.48795,96.1],[47.45382,8.55791,470],[47.59975,7.53157,314.7]]
Timestamps=[1623924000004404293,1623924000003618910,1623924000003779350]
Timestamps, Receivers = zip(*sorted(zip(Timestamps, Receivers)))
Timestamps=list(Timestamps)
Receivers=list(Receivers)
referencePx=Receivers[0]
referenceTimestamp=Timestamps[0]

Timestamps.remove(referenceTimestamp)
Receivers.remove(referencePx)

print(Timestamps,Receivers)

def rotateBack(X1,Y1):
    Y1R = X1 * math.sin(angleRotation) + Y1*math.cos(angleRotation)
    X1R = X1 * math.cos(angleRotation) - Y1*math.sin(angleRotation)
    return X1R,Y1R


lat=47.684189
lon=8.524155
alts=0

referenceCoord = GPS_utils()
referenceCoord.setENUorigin(referencePx[0],referencePx[1],referencePx[2])
theRealDeal=referenceCoord.geo2enu(lat,lon,alts*1000)
Xreal=theRealDeal[0,0]/1000
Yreal=theRealDeal[1,0]/1000

xi=0
yi=0
Zi=0-alts

xj=referenceCoord.geo2enu(Receivers[0][0],Receivers[0][1],Receivers[0][2])[0,0]/1000
yj=referenceCoord.geo2enu(Receivers[0][0],Receivers[0][1],Receivers[0][2])[1,0]/1000
Zj=referenceCoord.geo2enu(Receivers[0][0],Receivers[0][1],Receivers[0][2])[2,0]/1000-alts

xk=referenceCoord.geo2enu(Receivers[1][0],Receivers[1][1],Receivers[1][2])[0,0]/1000
yk=referenceCoord.geo2enu(Receivers[1][0],Receivers[1][1],Receivers[1][2])[1,0]/1000
Zk=referenceCoord.geo2enu(Receivers[1][0],Receivers[1][1],Receivers[1][2])[2,0]/1000-alts

angleRotation=(math.atan2(yk,xk))
#rotate the plane over the reference point
xkR = xk * math.cos(angleRotation) + yk * math.sin(angleRotation)
ykR = 0
xk=xkR
yk=ykR

xjR = xj * math.cos(angleRotation) + yj * math.sin(angleRotation)
yjR = -xj * math.sin(angleRotation) + yj * math.cos(angleRotation)
xj=xjR
yj=yjR

#now the axes have been rotated, start calculations

#define variables
c=300000000
Rij = c*((referenceTimestamp)-(Timestamps[0]))/1e9
Rik = c*((referenceTimestamp)-(Timestamps[1]))/1e9


DSj=xj**2+yj**2
DSk=xk**2
A=(Rij**2-DSj+Zi**2-(Zj**2))/2
B=(Rik**2-DSk+Zi**2-(Zk**2))/2


C=(Rij*xk-Rik*xj)/(Rik*yj)
D=(Rij*B - Rik*A) / (Rik*yj)

E=((xk/Rik)*(B/Rik) - C*D)/(C**2 - ((xk/Rik)**2) + 1) 
F = (D**2 + Zi**2 - (B/Rik)**2)/(C**2 - ((xk/Rik)**2) + 1)


#corner cases
###{...}


#normal case
X1 = E + (E**2 - F)**0.5
Y1 = C*X1 + D
X2 = E - (E**2 - F)**0.5
Y2 = C*X2 + D

X1,Y1=rotateBack(X1,Y1)
X2,Y2=rotateBack(X2,Y2)

print(X1,Y1)
print(X2,Y2)
print(Xreal,Yreal)