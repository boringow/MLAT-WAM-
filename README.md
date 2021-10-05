# MLAT-WAM-
Solve 2D Multilateration problem
    referenceRx=getReferenceRx(rec_times)
    referencePx=getPositionFromRx(referenceRx)
    referenceTime=min(rec_times)

    #We remove the reference from the list
    receivers.remove(referenceRx)
    rec_times.remove(referenceTime)

    if (len(receivers)==2):
        #order the receivers by time of arrival
        if (rec_times[0])>(rec_times[1]):
            for j in range(0,len(receivers)):
                RxLocations.append(getPositionFromRx(receivers[j-1]))
            rec_times[0],rec_times[1]=rec_times[1],rec_times[0]#swap the receiver time.
        else:
            for j in range(0,len(receivers)):
                RxLocations.append(getPositionFromRx(receivers[j]))# we dont need to swap receiver time as it is already in order of TOA


      def rotateBack(X1,Y1):
          Y1R = X1 * math.sin(angleRotation) + Y1*math.cos(angleRotation)
          X1R = X1 * math.cos(angleRotation) - Y1*math.sin(angleRotation)

          return X1R,Y1R

        referenceCoord = GPS_utils()
        referenceCoord.setENUorigin(referencePx[0],referencePx[1],referencePx[2])
        theRealDeal=referenceCoord.geo2enu(lat,lon,alts*1000)
        Xreal=theRealDeal[0,0]/1000
        Yreal=theRealDeal[1,0]/1000

        xi=0
        yi=0
        Zi=0-alts
        
        xj=referenceCoord.geo2enu(RxLocations[0][0],RxLocations[0][1],RxLocations[0][2])[0,0]/1000
        yj=referenceCoord.geo2enu(RxLocations[0][0],RxLocations[0][1],RxLocations[0][2])[1,0]/1000
        Zj=referenceCoord.geo2enu(RxLocations[0][0],RxLocations[0][1],RxLocations[0][2])[2,0]/1000-alts

        xk=referenceCoord.geo2enu(RxLocations[1][0],RxLocations[1][1],RxLocations[1][2])[0,0]/1000
        yk=referenceCoord.geo2enu(RxLocations[1][0],RxLocations[1][1],RxLocations[1][2])[1,0]/1000
        Zk=referenceCoord.geo2enu(RxLocations[1][0],RxLocations[1][1],RxLocations[1][2])[2,0]/1000-alts

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
        Rij = c*((referenceTime)-(rec_times[0]))/1e9
        Rik = c*((referenceTime)-(rec_times[1]))/1e9


        DSj=xj**2+yj**2
        DSk=xk**2
        A=(Rij**2-DSj+Zi**2-(Zj**2))/2
        B=(Rik**2-DSk+Zi**2-(Zk**2))/2


        C=(Rij*xk-Rik*xj)/(Rik*yj)
        D=(Rij*B - Rik*A) / (Rik*yj)
        
        E=((xk/Rik)*(B/Rik) - C*D)/(C**2 - ((xk/Rik)**2) + 1) 
        F = (D**2 + Zi**2 - (B/Rik)**2)/(C**2 - ((xk/Rik)**2) + 1)

        #corner cases
        if abs(yj)<3:#yj really small
            print("yj small")
            X1 = (Rik*A - Rij*B)/ (Rij*xk - Rik*xj)
            X2=X1

            K = (X1**2 + Zi**2)
            L = ((xj-X1)**2 + Zj**2)
            M = (Rij**2 + K - L) / (2*Rij)
            Y1 = (M**2 - K)**0.5
            Y2 = -(M**2 - K)**0.5

            X1,Y1=rotateBack(X1,Y1)
            X2,Y2=rotateBack(X2,Y2)
            continue

        if abs(Rik)<3:
            print("Rik small")
            X1 = -B / xk
            X2=X1
            if  abs(Rij)<3:
                print("Rij AND Rik small")
                Y1 = -(A + xj*X1) / yj
                Y2=Y1
                X1,Y1=rotateBack(X1,Y1)
                X2,Y2=rotateBack(X2,Y2)
                continue
            
            I = G*H / (1-H**2)
            J = ((B/xk)**2 + Zi**2 - G2) / (1-H**2)
            Y1 = I + (I**2 - J )**0.5
            Y2 = I - (I**2 - J )**0.5
            X1,Y1=rotateBack(X1,Y1)
            X2,Y2=rotateBack(X2,Y2)
            continue
        
        if (E**2 - F)<0:
            X1=E
            X2=X1
            Y1 = C*X1 + D
            Y2 = Y1
            X1,Y1=rotateBack(X1,Y1)
            X2,Y2=rotateBack(X2,Y2)
            continue

        #normal case
        X1 = E + (E**2 - F)**0.5
        Y1 = C*X1 + D
        X2 = E - (E**2 - F)**0.5
        Y2 = C*X2 + D

        X1,Y1=rotateBack(X1,Y1)
        X2,Y2=rotateBack(X2,Y2)
