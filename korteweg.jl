#= 
Korteweg de Vries Simulationsprogramm
=#
function PeriodicCond(pos)# Periodische Randbedingungen
  temp1 = u[pos,scal-3]
  temp2 = u[pos,scal-2]
  u[pos,1], u[pos,2],u[pos,scal-1], u[pos,scal] =
  temp1, temp2, u[pos,3], u[pos,4]
end
function Maximum(n)
  max = 0
  if n == 1
    for i in 3:(scal-2)
      if max < u[1,i]
        max = u[1,i]
      end
    end
  else
    for i in 3:(scal-2)
      if max < u[3,i]
        max = u[3,i]
      end
    end
  end
  return max
end
function Init()
  global file = open("P_K-data.txt","w")
  global P0alt =0.
  global P1alt =0.
  global P2alt =0.
  global printer =Int64[0:200:parse(Int64,ARGS[2]);]# Variable zur Steuerung der 
  #Ausgabe von Dateien, alle 200 Schritt erfolgt die Ausgabe
  printer[1]=1
  global schalter=parse(Int,ARGS[1])
  global N = parse(Float64,ARGS[5]) # N- Soliton Lösung
  global delta_x = parse(Float64,ARGS[3])
  global delta_t = 2*delta_x/(3*sqrt(3)*(-2*(N*(N+1))+1/(delta_x)^2))*0.7
  global grenze = parse(Float64,ARGS[4]) # Größe des Gitters
  global L = grenze*0.3
  global h0 = parse(Float64,ARGS[6]) # Tiefe des Wassers
  global l = Int64(2)
  global pil = Int(1)
  global x = Float64[-grenze:delta_x:3*grenze;] # x Koordinate
  global scal = length(x)+4# Indexgröße von u festlegen
  global u = (zeros(Float64,3,scal)) # dimensionslose Anhebung
  for i in 3:(scal-2) #Anfangsbedinung anlegen
    u[1,i] = ((N*(N+1))/(cosh(x[i-2]+grenze*(-parse(Float64,ARGS[8]))))^2)
  end
  println("reset")
  println("#$(parse(Int,ARGS[1]))-$(parse(Float64,ARGS[2]))-\\
  $(parse(Float64,ARGS[3]))-$(parse(Float64,ARGS[4]))-\\
  $(parse(Float64,ARGS[5]))-$(parse(Float64,ARGS[6]))-\\
  $(parse(Float64,ARGS[7]))-$(parse(Float64,ARGS[8]))")
  println("# Maxmimum $(Maximum(1))")
  println("set grid")
  println("set xrange [",-grenze,":",3*grenze,"]")
  println("set yrange [",-2,":",N*(N+5)/h0,"]")#4.5 , 10.5 , 18
  println("set xlabel 'x'")
  println("set ylabel 'u'")
end
function h(pos)
  if (schalter == 1) # h_1(x) Funktion
    if (x[pos] < 0)
      return 1
    elseif (0<=x[pos])&&(x[pos]<=L)
      return 1/2*(1+h0+(1-h0)*cos(pi*(x[pos])/L))
    else
      return h0
    end
  elseif (schalter == 2) # eine Testfunktion
    if (L<=x[pos])&&(x[pos]<=2*L)
      return 1/2*(1+h0+(1-h0)*cos(pi*(x[pos]+L)/L))
    elseif x[pos] > L
      return h0
    elseif (x[pos]<=-grenze+L)&&(x[pos]>=-grenze)
      return 1/2*(1+h0+(1-h0)*cos(pi*(x[pos]-grenze+L)/L))
    else
      return 1
    end
  elseif (schalter == 3) # h_2(x) Funktion
    return  (1-h0*exp(-x[pos]*x[pos]/(L*L)))
  else
    return 1
  end
end
function NumCalc(n)
        for i in 3:(scal-2) # Zabuski-Kruskal Schema
                if n == 2
                        u[2,i]=u[1,i] - delta_t*(((h(i-2))^(-7/4))*
                                (u[1,i-1]+u[1,i]+u[1,i+1])*
                                (u[1,i+1]-u[1,i-1])/delta_x+(h(i-2)^(1/2))*
                                (u[1,i+2]-2*u[1,i+1]+2*u[1,i-1]-u[1,i-2])/
                                (2*(delta_x)^3))
                        global l=2
                else
                        u[3,i]=u[1,i] - 2*delta_t*(((h(i-2))^(-7/4))*
                                (u[2,i-1]+u[2,i]+u[2,i+1])*
                                (u[2,i+1]-u[2,i-1])/delta_x+(h(i-2)^(1/2))*
                                (u[2,i+2]-2*u[2,i+1]+2*u[2,i-1]-u[2,i-2])/
                                (2*(delta_x)^3))
                        global l=3
                end
        end
end
function IntCon(pos,n)
  P0=Float64(0)
  for i in 3:(scal-2)
    P0=P0+(1/3*(u[pos,i-1]+u[pos,i]+u[pos,i+1]))*delta_x
  end
  P1=Float64(0)
  for i in 3:(scal-2)
    P1=P1+((1/3*(u[pos,i-1]+u[pos,i]+u[pos,i+1]))^2)*delta_x
  end
  P2=Float64(0)
  for i in 3:(scal-2)
    P2=P2+(2*(1/3*(u[pos,i-1]+u[pos,i]+u[pos,i+1]))^3-
      ((u[pos,i+1]-u[pos,i-1])/delta_x/2)^2)*delta_x
  end
  #println("set label 3 \"P0= ",
  # trunc(P0,digits=10)," \" at screen 0.6,0.8")# Anzeige des Integrals im Diagramm
  #println("set label 2 \"P1= ",
  # trunc(P1,digits=10)," \" at screen 0.6,0.75")# Anzeige des Integrals im Diagramm
  #println("set label 1 \"P2= ",
  # trunc(P2,digits=10)," \" at screen 0.6,0.7")# Anzeige des Integrals im Diagramm
  if n>1
    #write(file,"$n $(trunc(abs(Pk[1]-P0),digits=10)) $(trunc(abs(Pk[2]-P1),digits=10)) $(trunc(abs(Pk[3]-P2),digits=10)) \n")
    write(file,"$n $(abs(Pk[1]-P0)) $(abs(Pk[2]-P1)) $(abs(Pk[3]-P2)) \n")
    #write(file,"$n $(abs(Pk[1]-P0)) $(Pk[1]) \n")
  else
    global Pk=zeros(3)
    Pk[1]=P0
    Pk[2]=P1
    Pk[3]=P2
    write(file,"$n 0 0 0 \n")
    #write(file,"$n $(Pk[1]) $(Pk[2]) $(Pk[3]) \n")
  end
end
function Output(pos, tsim, treal, n) # Ausgabe ins Terminal
  if (pil <= length(printer)) & ( parse(Int,ARGS[7])> 0 )
    if (n == printer[pil]) 
      #println("set term qt;")
      println("set term png; set output sprintf('korteweg_frame$n-$(parse(Int,ARGS[1]))-$(parse(Float64,ARGS[3]))-$(parse(Float64,ARGS[4]))-$(parse(Float64,ARGS[5]))-$(parse(Float64,ARGS[6]))-$(parse(Float64,ARGS[8])).png');")
      global pil = pil+1
      println("set label 5 \"sim time= ",
      trunc(tsim,digits=10)," \" at screen 0.6,0.9")# Anzeige der Zeit im Diagramm
      println("set label 6 \"n= ",n," \" at screen 0.6,0.85")
      if (schalter == 1) # Anzeige der h_1,2 Funktionen im Diagramm
      println("plot '-' w l, [-$grenze:0] -1 linecolor 3 notitle, [0:$L] \\
      -(1+$(h0)+(1-$(h0))*cos(pi*x/$(L)))/2 linecolor 3 notitle, \\
      [$L:$(5*grenze)] -$(h0) linecolor 3 notitle")
      elseif (schalter == 2)
        println("plot '-' w l, [$L:$(2*L)] (1+$(h0)+(1-$(h0))*cos((pi*(x+$L))/$L))/2-1\\
        linecolor 3 notitle, [-$(grenze):$(-grenze+L)] (1+$(h0)+\\
        (1-$(h0))*cos(pi*(x-$(grenze+L))/$L))/2-1 linecolor 3 notitle, \\
        [$(-grenze+L):$L] 1-1 linecolor 3 notitle, [$(2*L):$grenze] $(h0)-1 linecolor 3 notitle") 
      elseif (schalter == 3)
        println("plot '-' w l, [-$grenze:$(3*grenze)] -(1-$(h0)*exp(-x*x/$(L*L))) linecolor 3 notitle")
      else
        println("plot '-' w l")
      end
      global y=zeros(scal-4)
      for i in 3:(scal-2) # Werteausgabe 
        println(x[i-2]," ",u[pos,i])
      end
      println("elapsed time: 0.672879991 seconds")
    end
  end
  if n == (parse(Int,ARGS[2])) # Maximum Anzeige und Berechnung
    println("# Maxmimum $(Maximum(n))")
  end
  #close("P_K.txt")
end
for n in 1:(parse(Int,ARGS[2]))
        if n == 1
                global treal=0.
                global tsim=0.
                Init()
                IntCon(n,n)
                PeriodicCond(n)
                tsim=delta_t
                Output(n, tsim, treal, n)
        else
                NumCalc(n)
                PeriodicCond(l)
                if n>2 #Vertauschung der Reihen in dem Array der Anhebung
                        for k in 1:scal
                                temp3=u[2,k]
                                u[2,k]=u[3,k]
                                u[1,k]=temp3
                        end
                end
                tsim=tsim+delta_t
                #println("elapsed time: 0.672879991 seconds")
                IntCon(l,n)
                Output(l, tsim, treal, n)
        end
end

