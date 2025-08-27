#= Programm für das Computational Physics Projekt 
   Korteweg de Vries Equation
=#

# Problem 1

global P0alt =0.
global P1alt =0.
global P2alt =0.

function PeriodicCond(pos)# Periodische Randbedingungen
    # Functions return the value of their last statement
    temp1 = u[pos,scal-3]
    temp2 = u[pos,scal-2]
    u[pos,1], u[pos,2],u[pos,scal-1], u[pos,scal] = temp1, temp2, u[pos,3], u[pos,4]
end

function Init()
    N = 2 # N- soliton Lösung
    global delta_x = 0.05
    global delta_t = 2*delta_x/(3*sqrt(3)*(-2*(N*(N+1))+1/(delta_x)^2))*0.7
    #global delta_t =0.0001
    println("# ", delta_t)
    global grenze = 5
    global L = grenze/6
    global h0 = N*0.4
    global l = 2
    global x = [-grenze:delta_x:grenze;] # x Koordinate
    global scal = length(x)+4
    #println(scal)
    global u = zeros(3,scal)
    #println(x)
    println("# ", 1)
    for i in 3:(scal-2)
        u[1,i] = ((N*(N+1))/(cosh(x[i-2]+1*grenze*0.3))^2)
        #println("#", x[i-2]," ",u[1,i])
    end
    println("set xrange [",-grenze,":",grenze,"]")
    println("set yrange [",-1,":",N*(N+1)+N*2,"]")
end

function h(pos)
    if (L<=x[pos])&&(x[pos]<=2*L)
        return 1./2.*(1+h0+(1-h0)*cos(pi*(x[pos]+L)/L))
    elseif x[pos] > L
        return 1
    elseif (x[pos]<=-grenze+L)&&(x[pos]>=-grenze)
        return 1./2.*(1+h0+(1-h0)*cos(pi*(x[pos]-grenze+L)/L))
    #elseif x[pos]<-2*L
    #    return 1
    else
        return h0
    end
end

function NumCalc(n)
        for i in 3:(scal-2)
            if n == 2
                u[2,i]=u[1,i] - delta_t*(((h(i-2))^(-7/4))*(u[1,i-1]+u[1,i]+u[1,i+1])*(u[1,i+1]-u[1,i-1])/delta_x+(h(i-2)^(1/2))*(u[1,i+2]-2*u[1,i+1]+2*u[1,i-1]-u[1,i-2])/(2*(delta_x)^3))
                l = 2
            else
                u[3,i]=u[1,i] - 2*delta_t*(((h(i-2))^(-7/4))*(u[2,i-1]+u[2,i]+u[2,i+1])*(u[2,i+1]-u[2,i-1])/delta_x+(h(i-2)^(1/2))*(u[2,i+2]-2*u[2,i+1]+2*u[2,i-1]-u[2,i-2])/(2*(delta_x)^3))
                l = 3
            end
        end
end

function IntCon(pos)
    P0=0.
    for i in 3:(scal-2)
        P0=P0+u[pos,i]*delta_x
    end
    P1=0.
    for i in 3:(scal-2)
        P1=P1+(u[pos,i]^2)*delta_x
    end
    P2=0.
    for i in 3:(scal-2)
        P2=P2+(2*(u[pos,i])^3-((u[pos,i+1]-u[pos,i-1])/delta_x)^2)*delta_x
    end
    println("set label 1 \"P0= ",P0," \" at screen 0.6,0.8")# Anzeige des Integrals im Diagramm
    println("set label 2 \"P1= ",P1," \" at screen 0.6,0.75")# Anzeige des Integrals im Diagramm
    println("set label 3 \"P2= ",P2," \" at screen 0.6,0.7")# Anzeige des Integrals im Diagramm
    #println("#P0 ", P0-P0alt, " P1=", P1-P1alt," P2=", P2-P2alt)
    #P0alt=P0
    #P1alt=P1
    #P2alt=P2
end

function Output(pos, tsim, treal) # Ausgabe ins Terminal
    println("set label 4 \"sim time= ",tsim," \" at screen 0.6,0.9")# Anzeige der Zeit im Diagramm
    println("set label 5 \"real time= ",treal," \" at screen 0.6,0.85")# Anzeige der Zeit im Diagramm
    println("plot '-' w l")
    for i in 3:(scal-2)
        println(x[i-2]," ",u[pos,i])
    end
     
#= 
    open("t.txt", "w") do f
        for i in 3:(scal-2)
            x[i-2], u[3,i]
            write(f,"$x,  $u \n")
        end
    end
=#

end

T = 50000
treal=0.
tsim=0.
for n in 1:T
    tic()# start timer und Anker für den nächste Messung; Zeit wird zwischen tic() und toc() gemessen und bei toc() ausgegeben
    if n == 1
        Init()
        IntCon(n)
        PeriodicCond(n)
        tsim=delta_t
        Output(n, tsim, treal)
    else
        NumCalc(n)
        PeriodicCond(l)
        if n>2
            for k in 1:scal
                temp1=u[2,k]
                u[2,k]=u[3,k]
                u[1,k]=temp1
            end
        end
        tsim=tsim+delta_t
        treal=treal+toc()#Ende der Messung übertragen an den Zähler zur Ausgabe
        #sleep(0.05)
        IntCon(l)
        Output(l, tsim, treal)
    end 
end
