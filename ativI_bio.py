import matplotlib.pyplot as plt
import numpy as np 
from math import e,pi,sqrt
from matplotlib.ticker import NullFormatter


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

def histogramas():

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

    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
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

    if aminoacido == "L" or aminoacido=="l":

        plt.title('Leucina')
        plt.xlabel('Phi(φ)')
        plt.ylabel('Frequência Absoluta')
        plt.hist(PhiLeu_33, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')

        plt.show()

        plt.title('Leucina')
        plt.xlabel('Psi(ψ)')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiLeu_33, bin, rwidth=2, color='red', alpha=1,edgecolor='black')
        
        plt.show()

    elif aminoacido == "T" or aminoacido=="t":
        plt.title('Theta')
        plt.xlabel('Theta(Θ)')
        plt.ylabel('Frequência Absoluta')
        
        plt.hist(angulo, bin, rwidth=2, color='blue', alpha=0.7,edgecolor='black')

        plt.show()

    else:
        plt.title('Ácido glutâmico')
        plt.xlabel('Phi(φ)')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PhiGlu_10, bin, rwidth=0.9, color='green', alpha=0.7,edgecolor='black')

        plt.show()

        plt.title('Ácido glutâmico')
        plt.xlabel('Psi(ψ)')
        plt.ylabel('Frequência Absoluta')

        plt.hist(PsiGlu_10, bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

        plt.show()

def funcaoProb():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []


    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
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

    # Phi e Psi Leu
    plt.title('Leucina')
    plt.xlabel('Phi(ψ)')
    plt.ylabel('Probabilidade')

    x = np.arange(-180,180,1)
    y = 1/(dp3*sqrt(2*pi))*np.exp(-1/2*((x-E3)/dp3)**2)
    plt.plot(x,y)

    plt.show()
    plt.title('Leucina')
    plt.xlabel('Psi(φ)')
    plt.ylabel('Probabilidade')

    x = np.arange(-180,180,1)
    y = 1/(dp4*sqrt(2*pi))*np.exp(-1/2*((x-E4)/dp4)**2)
    plt.plot(x,y)


    plt.show()

    plt.title('Theta')
    plt.xlabel('Theta(Θ)')
    plt.ylabel('Probabilidade')

    x = np.arange(-180,180,1)
    y = 1/(dp0*sqrt(2*pi))*np.exp((-1/2)*((x-E0)/dp0)**2)

    plt.plot(x,y)

    plt.show()

    # Phi e Psi Glu
    plt.title('Ácido Glutamínico')
    plt.xlabel('Phi(ψ)')
    plt.ylabel('Probabilidade')

    x = np.arange(-180,180,1)
    y = 1/(dp1*sqrt(2*pi))*np.exp(-1/2*((x-E1)/dp1)**2)
    plt.plot(x,y)

    plt.show()

    plt.title('Ácido Glutamínico')
    plt.xlabel('Psi(φ)')
    plt.ylabel('Probabilidade')

    x = np.arange(-180,180,1)
    y = 1/(dp2*sqrt(2*pi))*np.exp(-1/2*((x-E2)/dp2)**2)
    plt.plot(x,y)

    plt.show()

    arq.close()

def hist_2750_psixtempo():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
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
    plt.xlabel('Ângulo Psi')
    plt.ylabel('Frequência Absoluta')

    plt.hist(PsiGlu_10[2750:], bin, rwidth=0.9, color='purple', alpha=0.6,edgecolor='black')

    plt.show()


    plt.ylabel('Ângulo')
    plt.xlabel('Tempo (ps)')
    plt.plot(tempo,PsiGlu_10)
    plt.plot(tempo,angulo)

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

    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
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
    print("Variância Theta: ",dp0)
    print("Variância Phi Glutamato: ",dp1)
    print("Variância Psi Glutamato: ",dp2)
    print("Variância Phi Leucina: ",dp3)
    print("Variância Psi Leucina: ",dp4)
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

    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
            matriz.append(list(map(float,linha.strip("\n").split(","))))


    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    plt.title('Ramachandran Leu')
    plt.xlabel('Phi(ψ)')
    plt.ylabel('Psi(φ)')
    plt.hist2d(PhiLeu_33,PsiLeu_33,bins=(150,150), range=[(-180,180), (-180,180)])
    plt.axhline(y=0,linestyle="solid")
    plt.axvline(x=0,linestyle="solid")
    plt.show()

    plt.title('Ramachandran Glu')
    plt.xlabel('Phi(ψ)')
    plt.ylabel('Psi(φ)')
    plt.hist2d(PhiGlu_10,PsiGlu_10,bins=(150,150), range=[(-180,180), (-180,180)])
    plt.axhline(y=0,linestyle="solid")
    plt.axvline(x=0,linestyle="solid")
    plt.show()

def projecao():

    arq = open("trajetorias.csv","r")

    matriz = []

    tempo = []
    angulo = []
    PhiGlu_10 = []
    PsiGlu_10 = []
    PhiLeu_33 = []
    PsiLeu_33 = []

    bin = []

    i = 1
    for linha in arq:
        if i ==1:
            i+=1
            pass
        else:
            matriz.append(list(map(float,linha.strip("\n").split(","))))

    for linha in matriz:
        tempo.append(linha[0])
        angulo.append(linha[1])
        PhiGlu_10.append(linha[2])
        PsiGlu_10.append(linha[3])
        PhiLeu_33.append(linha[4])
        PsiLeu_33.append(linha[5])

    mat = [PhiGlu_10,PsiGlu_10,PhiLeu_33,PsiLeu_33]
    for i in range(2):

        x = mat[i%2*2]
        y = mat[i%2*2+1]
        nullfmt = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure
        plt.figure(1, figsize=(8, 8))

        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        # the scatter plot:
        axScatter.scatter(x, y,0.5)

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
        lim = (int(xymax/binwidth) + 1) * binwidth

        axScatter.set_xlim((-lim, lim))
        axScatter.set_ylim((-lim, lim))

        bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(x, bins=bins)
        axHisty.hist(y, bins=bins, orientation='horizontal')

        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())

        plt.show()
