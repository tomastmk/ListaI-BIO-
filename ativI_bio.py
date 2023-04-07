import matplotlib.pyplot as plt
import numpy as np 
from math import e,pi,sqrt

def esperanca(lista,n):
    E = 0
    for i in lista:
        E +=i
    E /= n
    return E

def dpAmostral(lista,esperanca,n):
    soma = 0
    for i in lista:
        soma = soma +(i-esperanca)**2
    dp = (soma/(n-1))**(1/2)
    return dp

def aI():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    aminoacido = input("Proteína (G ou L) ou Theta (T):  ")

    for linha in arq:
        matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    E0 = esperanca(angulo,5500)
    E1 = esperanca(PhiGlu_10,5500)
    E2 = esperanca(PsiGlu_10,5500)
    E3 = esperanca(PhiLeu_33,5500)
    E4 = esperanca(PsiLeu_33,5500)
    dp0 = dpAmostral(angulo,E0,5500)
    dp1 = dpAmostral(PhiGlu_10,E1,5500)
    dp2 = dpAmostral(PsiGlu_10,E2,5500)
    dp3 = dpAmostral(PhiLeu_33,E3,5500)
    dp4 = dpAmostral(PsiLeu_33,E4,5500)

    for i in range(40):
        bin.append(-180+(i+1)*4.5)
    for i in range(40):
        bin.append((i+1)*4.5)

    if aminoacido == "L" or aminoacido=="l":

        plt.title('Leucina')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')
        plt.hist(PhiLeu_33, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')

        plt.show()

        plt.title('Leucina')
        plt.xlabel('Ângulo Psi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiLeu_33, bin, rwidth=2, color='red', alpha=1,edgecolor='black')
        
        plt.show()

    elif aminoacido == "T" or aminoacido=="t":
        plt.title('Theta')
        plt.xlabel('Ângulo Theta')
        plt.ylabel('Frequência Absoluta')
        
        plt.hist(angulo, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')

        plt.show()

    else:
        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PhiGlu_10, bin, rwidth=0.9, color='green', alpha=0.7,edgecolor='black')

        plt.show()

        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiGlu_10, bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

        plt.show()

def aII():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []


    for linha in arq:
        matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    E0 = esperanca(angulo,5500)
    E1 = esperanca(PhiGlu_10,5500)
    E2 = esperanca(PsiGlu_10,5500)
    E3 = esperanca(PhiLeu_33,5500)
    E4 = esperanca(PsiLeu_33,5500)
    dp0 = dpAmostral(angulo,E0,5500)
    dp1 = dpAmostral(PhiGlu_10,E1,5500)
    dp2 = dpAmostral(PsiGlu_10,E2,5500)
    dp3 = dpAmostral(PhiLeu_33,E3,5500)
    dp4 = dpAmostral(PsiLeu_33,E4,5500)

    print()
    print("E theta: ",E0)
    print("E Phi Glutamato: ",E1)
    print("E Psi Glutamato: ",E2)
    print("E Phi Leucina: ",E3)
    print("E Psi Leucina: ",E4)
    print()
    print()
    print("DP Theta: ",dp0)
    print("DP Phi Glutamato: ",dp1)
    print("DP Psi Glutamato: ",dp2)
    print("DP Phi Leucina: ",dp3)
    print("DP Psi Leucina: ",dp4)
    print()
    print("%s%f%s%f%s%f%s" %("1/(",dp0,"*sqrt(2*pi))*np.exp(-1/2*((x-",E0,")/",dp0,")**2)"))
    print("%s%f%s%f%s%f%s" %("1/(",dp1,"*sqrt(2*pi))*np.exp(-1/2*((x-",E1,")/",dp1,")**2)"))
    print("%s%f%s%f%s%f%s" %("1/(",dp2,"*sqrt(2*pi))*np.exp(-1/2*((x-",E2,")/",dp2,")**2)"))
    print("%s%f%s%f%s%f%s" %("1/(",dp3,"*sqrt(2*pi))*np.exp(-1/2*((x-",E3,")/",dp3,")**2)"))
    print("%s%f%s%f%s%f%s" %("1/(",dp4,"*sqrt(2*pi))*np.exp(-1/2*((x-",E4,")/",dp4,")**2)"))
    print()


    for i in range(40):
        bin.append(-180+(i+1)*4.5)
    for i in range(40):
        bin.append((i+1)*4.5)

    aminoacido = input("Proteína (G ou L) ou Theta (T):  ")

    if aminoacido == "L" or aminoacido=="l":

        plt.title('Leucina')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')
        plt.hist(PhiLeu_33, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')
        
        x = np.arange(-180,180,1)
        y = 1/(dp3*sqrt(2*pi))*np.exp(-1/2*((x-E3)/dp3)**2)*78000
        plt.plot(x,y)
        plt.show()


        plt.title('Leucina')
        plt.xlabel('Ângulo Psi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiLeu_33, bin, rwidth=2, color='red', alpha=1,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp4*sqrt(2*pi))*np.exp(-1/2*((x-E4)/dp4)**2)*40725
        plt.plot(x,y)


        plt.show()

    elif aminoacido == "T" or aminoacido=="t":
        plt.title('Theta')
        plt.xlabel('Ângulo Theta')
        plt.ylabel('Frequência Absoluta')

        plt.hist(angulo, bin, rwidth=2, color='orange', alpha=0.7,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp0*sqrt(2*pi))*np.exp((-1/2)*((x-E0)/dp0)**2)*22807

        plt.plot(x,y)

        plt.show()

    else:
        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PhiGlu_10, bin, rwidth=0.9, color='green', alpha=0.7,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp1*sqrt(2*pi))*np.exp(-1/2*((x-E1)/dp1)**2)*27733
        plt.plot(x,y)

        plt.show()
        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiGlu_10, bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp2*sqrt(2*pi))*np.exp(-1/2*((x-E2)/dp2)**2)*51000
        plt.plot(x,y)

        plt.show()
        arq.close()

def aIII():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    for linha in arq:
        matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])


    for i in range(40):
        bin.append(-180+(i+1)*4.5)
    for i in range(40):
        bin.append((i+1)*4.5)

    plt.title('Theta')
    plt.xlabel('Ângulo Theta')
    plt.ylabel('Frequência Absoluta')
    
    plt.hist(angulo[2750:], bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')

    plt.show()


    plt.title('Glutamato')
    plt.xlabel('Ângulo Phi')
    plt.ylabel('Frequência Absoluta')

    plt.hist(PhiGlu_10[2750:], bin, rwidth=0.9, color='green', alpha=0.7,edgecolor='black')

    plt.show()

    plt.title('Glutamato')
    plt.xlabel('Ângulo Phi')
    plt.ylabel('Frequência Absoluta')

    plt.hist(PsiGlu_10[2750:], bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

    plt.show()


    plt.plot(tempo,angulo)
    plt.plot(tempo,PhiGlu_10)
    plt.plot(tempo,PsiGlu_10)

    plt.show()

def funcao():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    aminoacido = input("Proteína (G ou L) ou Theta (T):  ")

    for linha in arq:
        matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    E0 = esperanca(angulo,5500)
    E1 = esperanca(PhiGlu_10,5500)
    E2 = esperanca(PsiGlu_10,5500)
    E3 = esperanca(PhiLeu_33,5500)
    E4 = esperanca(PsiLeu_33,5500)
    dp0 = dpAmostral(angulo,E0,5500)
    dp1 = dpAmostral(PhiGlu_10,E1,5500)
    dp2 = dpAmostral(PsiGlu_10,E2,5500)
    dp3 = dpAmostral(PhiLeu_33,E3,5500)
    dp4 = dpAmostral(PsiLeu_33,E4,5500)

    print()
    print("E theta: ",E0)
    print("E Phi Glutamato: ",E1)
    print("E Psi Glutamato: ",E2)
    print("E Phi Leucina: ",E3)
    print("E Psi Leucina: ",E4)
    print()
    print()
    print("DP Theta: ",dp0)
    print("DP Phi Glutamato: ",dp1)
    print("DP Psi Glutamato: ",dp2)
    print("DP Phi Leucina: ",dp3)
    print("DP Psi Leucina: ",dp4)
    print()

    for i in range(40):
        bin.append(-180+(i+1)*4.5)
    for i in range(40):
        bin.append((i+1)*4.5)

    if aminoacido == "L" or aminoacido=="l":

        plt.title('Leucina')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')
        plt.hist(PhiLeu_33, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')
        
        x = np.arange(-180,180,1)
        y = 1/(dp3*sqrt(2*pi))*np.exp(-1/2*((x-E3)/dp3)**2)
        plt.plot(x,y)
        plt.show()


        plt.title('Leucina')
        plt.xlabel('Ângulo Psi')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiLeu_33, bin, rwidth=2, color='red', alpha=1,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp4*sqrt(2*pi))*np.exp(-1/2*((x-E4)/dp4)**2)
        plt.plot(x,y)


        plt.show()

    elif aminoacido == "T" or aminoacido=="t":
        plt.title('Theta')
        plt.xlabel('Ângulo Theta')
        plt.ylabel('Frequência Absoluta')
        plt.hist(angulo, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')
        x = np.arange(-180,180,1)
        y = 1/(dp0*sqrt(2*pi))*np.exp((-1/2)*((x-E0)/dp0)**2)
        plt.plot(x,y)
        plt.show()

    else:
        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        #plt.hist(PhiGlu_10, bin, rwidth=0.9, color='green', alpha=0.7,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp1*sqrt(2*pi))*np.exp(-1/2*((x-E1)/dp1)**2)
        plt.plot(x,y)

        plt.show()
        plt.title('Glutamato')
        plt.xlabel('Ângulo Phi')
        plt.ylabel('Frequência Absoluta')

        #plt.hist(PsiGlu_10, bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

        x = np.arange(-180,180,1)
        y = 1/(dp2*sqrt(2*pi))*np.exp(-1/2*((x-E2)/dp2)**2)
        plt.plot(x,y)

        plt.show()
        arq.close()

def ramachandran():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    for linha in arq:
        matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    plt.title('Glutamina')
    plt.xlabel('Ângulo Phi')
    plt.ylabel('Ângulo PSi')
    plt.scatter(PhiGlu_10,PsiGlu_10,0.5)
    
    plt.show()

    plt.title('Leucina')
    plt.xlabel('Ângulo Phi')
    plt.ylabel('Ângulo PSi')
    plt.scatter(PhiLeu_33,PsiLeu_33,0.5)
    
    plt.show()


ramachandran()

