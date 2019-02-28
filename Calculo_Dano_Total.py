# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:59:09 2019

@author: luize
"""

import numpy as np
#import pandas as pd

#Propriedades do material

Tensao_ultima = 25.76
a = 32.24
b = -0.068

#Definicao de ciclos para vida infinita

N_infinito = 1e8

#Definicao do fator de correcao de tensao alternada por Walker

gama = 0.4

#Fatores que alteram a vida com relacao a fadiga

Fator_de_Confiabilidade = 0.814

Fatores = Fator_de_Confiabilidade

#Define dados para expectativa de vida

Meta_de_vida = 250000
km_por_ciclo = 25
Expectativa_vida = Meta_de_vida / km_por_ciclo

#Funcao para importar valores pelo Python

def Extrair_Dados(Nome_Arquivo):
    csv = np.genfromtxt ('Dados\Fadiga\Proposta\%s'%(Nome_Arquivo)+'.csv', delimiter="\t",skip_header=1)
    return csv

#Funcao para encontrar dimensao de array e retorna valor inteiro
    
def Tamanho_Array(array):
    array_shape = array.shape
    Dim_array = int(array_shape[0])
    return Dim_array

#Calculo da maxima principal para cada no e seu estado de tensao em cada momento
    
def Calculo_Maxima(Vetor_Estado_Tensao,Numero_de_Nos):
    
    Vetor_Max_Principal = np.zeros((Numero_de_Nos,1))
    
    for i in range(0,Numero_de_Nos,1):
        
        Estado_Tensao = np.array(Vetor_Estado_Tensao[i,3:9],ndmin = 2)
        Tensor = np.zeros((3,3))
        
        for m in range (0,3,1):
            
            for n in range (0,3,1):
                
                if m == n:
                    if m == 0:
                        Tensor[m,n] = Estado_Tensao[0,m]
                    if m == 1:
                        Tensor[m,n] = Estado_Tensao[0,m]
                    if m == 2:
                        Tensor[m,n] = Estado_Tensao[0,m]
                if n > m:
                    if m < 1:
                        if n == 1:
                            Tensor[m,n] = Estado_Tensao[0,n+2]
                            Tensor[n,m] = Estado_Tensao[0,n+2]
                        else:
                            Tensor[m,n] = Estado_Tensao[0,n+3]
                            Tensor[n,m] = Estado_Tensao[0,n+3]
                    if m == 1:
                        Tensor[m,n] = Estado_Tensao[0,n+2]
                        Tensor[n,m] = Estado_Tensao[0,n+2]
                    
                    
        Max_Principais, Angulos = np.linalg.eigh(Tensor)
        
        MaxPrincipal1 = np.max(Max_Principais)
        
        Vetor_Max_Principal[i] = MaxPrincipal1

    return Vetor_Max_Principal

#Funcao para fazer contagem de ciclos em cada no e retornar dano, vida, fator de segurança e tensao equivalente

def Calcular_Fadiga(No):
    
    Ciclo_no = No

    #Ordena Ciclo de forma que tensao maxima fique nas extremidades
    
    D_Ciclo = Tamanho_Array(Ciclo_no)
    
    #Obtem posicao de maxima tensao no ciclo
    Posicao_max = np.argmax(Ciclo_no)
    
    #Dividi o ciclo em antes e depois do máximo
    
    Ciclo_no_antes = Ciclo_no[0:Posicao_max]
    Ciclo_no_depois = Ciclo_no[Posicao_max:D_Ciclo]
    Ciclo_Final = np.array(Ciclo_no[Posicao_max],ndmin=1)
    
    #Cria ciclo ordenado para calcular o Rainflow
    
    Ciclo_no_Ordenado = np.concatenate((Ciclo_no_depois,Ciclo_no_antes,Ciclo_Final))
   # D_Ciclo_no_Ordenado = Tamanho_Array(Ciclo_no_Ordenado)
    
    #Filtrar Ciclo_no_ordenado para manter apenas picos e vales
    
    Ciclo_no_Ordenado_pv = np.copy(Ciclo_no_Ordenado)
    
    D_Ciclo_no_Ordenado_pv = Tamanho_Array(Ciclo_no_Ordenado_pv)
    
   # D_Ciclo_no_Ordenado_pv_loop = D_Ciclo_no_Ordenado_pv
    
    #Cria vetor para identificar os carregamentos intermediarios que serao deletados
    Posicao_Deletar = np.zeros((D_Ciclo_no_Ordenado_pv,1),dtype=int)
    
    for i in range(1,D_Ciclo_no_Ordenado_pv-1,1):
        
        anterior = Ciclo_no_Ordenado_pv[i]-Ciclo_no_Ordenado_pv[i-1]
        proximo = Ciclo_no_Ordenado_pv[i+1]-Ciclo_no_Ordenado_pv[i]
        
        if anterior == 0 or proximo == 0 or np.sign(anterior) == np.sign(proximo):
            Ciclo_no_Ordenado_pv[i] = Ciclo_no_Ordenado_pv[i-1]
            Posicao_Deletar[i] = i
            
    Posicao_Deletar = np.array(np.nonzero(Posicao_Deletar),dtype=int, ndmin=0)
    
    Ciclo_no_Ordenado_pv = np.delete(Ciclo_no_Ordenado_pv, Posicao_Deletar[0])
        
        
    #Contar os ciclos e montar a matriz com as tensoes
    
    Ciclo_rainflow = np.copy(Ciclo_no_Ordenado_pv)
    
    D_rainflow = Tamanho_Array(Ciclo_rainflow)
    
    Tamanho_Matriz_Ciclos = int((D_rainflow-1)/2)
    
    Matriz_Ciclos = np.zeros((Tamanho_Matriz_Ciclos,8))
    Conta_Ciclos = 0
    
    #Vasculha os ciclos
    
    while D_rainflow>3:
        for i in range (1,D_rainflow-2,1):
            
            Ciclo_Externo = [Ciclo_rainflow[i-1],Ciclo_rainflow[i+2]]
            Ciclo_Interno = [Ciclo_rainflow[i],Ciclo_rainflow[i+1]]
            Delta_Externo = np.abs(Ciclo_Externo[0]-Ciclo_Externo[1])
            Delta_Interno = np.abs(Ciclo_Interno[0]-Ciclo_Interno[1])
            
            if Delta_Interno <= Delta_Externo :
                
                Minimo_Externo = Ciclo_Externo[np.argmin(Ciclo_Externo)]
                Maximo_Externo = Ciclo_Externo[np.argmax(Ciclo_Externo)]
                Minimo_Interno = Ciclo_Interno[np.argmin(Ciclo_Interno)]
                Maximo_Interno = Ciclo_Interno[np.argmax(Ciclo_Interno)]
                
                if Minimo_Interno >= Minimo_Externo and Maximo_Interno <= Maximo_Externo:
                    
                    Matriz_Ciclos[Conta_Ciclos,0] = Minimo_Interno
                    Matriz_Ciclos[Conta_Ciclos,1] = Maximo_Interno
                    
                    Conta_Ciclos = Conta_Ciclos + 1
                    
                    Ciclo_rainflow = np.delete(Ciclo_rainflow,(i,i+1))
                    break
                
        D_rainflow = Tamanho_Array(Ciclo_rainflow)
    
    #Ultimo ciclo
        Matriz_Ciclos[Conta_Ciclos,0] = Ciclo_rainflow[np.argmin(Ciclo_rainflow)]
        Matriz_Ciclos[Conta_Ciclos,1] = Ciclo_rainflow[np.argmax(Ciclo_rainflow)]
        
    #Calcula os valores para fadiga
    for i in range(0,Tamanho_Matriz_Ciclos,1):
        
        #Cálculo Tensão Média
        Matriz_Ciclos[i,2] = (Matriz_Ciclos[i,0]+Matriz_Ciclos[i,1])/2
        
        #Cálculo Tensão Alternada
        Matriz_Ciclos[i,3] = (Matriz_Ciclos[i,1]-Matriz_Ciclos[i,0])/2
        
        #Cálculo da Tensão Alternada corrigida
        
        
        #Cálculo Tensão ALternada corrigida pela Tensão Média Goodman
        
        #Matriz_Ciclos[i,4] = (Matriz_Ciclos[i,3]*Tensao_ultima/(Tensao_ultima-Matriz_Ciclos[i,2])) / Fatores
        
        #Cálculo Tensão ALternada corrigida pela Tensão Média Goodman para R=0
        
        #Matriz_Ciclos[i,4] = (Matriz_Ciclos[i,3]*Tensao_ultima/( (Matriz_Ciclos[i,3] - Matriz_Ciclos[i,2]) + Tensao_ultima)) / Fatores
        
        #Cálculo Tensão Alternada corrigida por Gerber
        
        #Matriz_Ciclos[i,4] = (Matriz_Ciclos[i,3]/(1-(Matriz_Ciclos[i,2]/Tensao_ultima)**2)) / Fatores
        
        #Cálculo Tensão Alternada corrigida por Walker
        
        Matriz_Ciclos[i,4] = ((Matriz_Ciclos[i,3]+Matriz_Ciclos[i,2])**(1-gama)*Matriz_Ciclos[i,3]**gama) / Fatores
        
        #Sem correção
        
        #Matriz_Ciclos[i,4] = Matriz_Ciclos[i,3] / Fatores
        
        #Cálculo Número de Ciclos
        Numero_Ciclos = (Matriz_Ciclos[i,4] / a) ** (1/b)
        if Numero_Ciclos < N_infinito:
            Matriz_Ciclos[i,5] = Numero_Ciclos
            Matriz_Ciclos[i,6] = 1 / Numero_Ciclos
        else:
            Matriz_Ciclos[i,5] = N_infinito
            Matriz_Ciclos[i,6] = 0
    
    #Regra de Palmgren Miner somar todos os danos
            
    #Soma os danos de todos os ciclos e multiplica pelo numero de repeticoes
            
    N_repeticoes=1
     
    Dano_Acumulado_no = 0
    
    for i in range(0,Tamanho_Matriz_Ciclos,1):
        Dano_Acumulado_no = Dano_Acumulado_no + Matriz_Ciclos[i,6]
           
    Dano_Acumulado_no = Dano_Acumulado_no*N_repeticoes
    
    #Verifica com base na quantidade desses ciclos quanto a estrutura suporta
    
    if Dano_Acumulado_no == 0:
        Vida_no = N_infinito
        Tensao_alternada_eq = 0
        Coef_Seg = 5
    else:
        if 1/Dano_Acumulado_no < N_infinito:
            
            Vida_no = 1/Dano_Acumulado_no
            
            #Calculada Tensao Equivalente
            
            Tensao_alternada_eq = a*Vida_no**b
            
            #Calcula o coeficiente de segurança
            
            if (Vida_no / Expectativa_vida) < 5:
                Coef_Seg = (Vida_no / Expectativa_vida)
            else:
                Coef_Seg = 5
            
        else:
            Vida_no = N_infinito
            Tensao_alternada_eq = 0
            Coef_Seg = 5
            
    return (Dano_Acumulado_no, Vida_no, Tensao_alternada_eq, Coef_Seg)
#Criando matrizes de tensão para cada etapa

Desmontado_Equilibrado = Extrair_Dados('Ciclos\Desmontado_Equilibrado')

D = Tamanho_Array(Desmontado_Equilibrado)

Montagem = Extrair_Dados('Ciclos\Montagem_Python_Ordenada')

Montado_Equlibrado = Extrair_Dados('Ciclos\Montado_Equilibrado')

Abrir_Quebrasol =  Extrair_Dados('Ciclos\Abrir_Quebrasol')

Fechar_Quebrasol =  Extrair_Dados('Ciclos\Fechar_Quebrasol')

Desmontagem =  Extrair_Dados('Ciclos\Desmontagem')

#Calcular Maxima Tensao Principal para cada array

Desmontado_Equilibrado_Max = np.zeros((D,1))
Desmontado_Equilibrado_Max = Calculo_Maxima(Desmontado_Equilibrado,D)

Montagem_Max = np.zeros((D,1))
Montagem_Max = Calculo_Maxima(Montagem,D)

Montado_Equilibrado_Max = np.zeros((D,1))
Montado_Equilibrado_Max = Calculo_Maxima(Montado_Equlibrado,D)

Abrir_Quebrasol_Max = np.zeros((D,1))
Abrir_Quebrasol_Max = Calculo_Maxima(Abrir_Quebrasol,D)

Fechar_Quebrasol_Max = np.zeros((D,1))
Fechar_Quebrasol_Max = Calculo_Maxima(Fechar_Quebrasol,D)

Desmontagem_Max = np.zeros((D,1))
Desmontagem_Max = Calculo_Maxima(Desmontagem,D)

#Montar Ciclo de fadiga

Montar = 1
Fechar_e_Abrir = 3
Desmontar = 1

Ciclo_Fadiga = Desmontado_Equilibrado_Max

for i in range(0,Montar,1):
    Ciclo_Fadiga = np.concatenate((Ciclo_Fadiga, Montagem_Max, Montado_Equilibrado_Max),axis =1)

for i in range(0,Fechar_e_Abrir,1):
    Ciclo_Fadiga = np.concatenate((Ciclo_Fadiga, Fechar_Quebrasol_Max, Montado_Equilibrado_Max, Abrir_Quebrasol_Max, Montado_Equilibrado_Max),axis =1)

for i in range(0,Desmontar,1):
    Ciclo_Fadiga = np.concatenate((Ciclo_Fadiga, Desmontagem_Max, Desmontado_Equilibrado_Max),axis =1)
    
#Obter dados fadiga cada no
    
Resultados_Fadiga_Ciclos = np.zeros((D,4))
    
for i in range(0,D,1):
    Resultados_Fadiga_Ciclos[i] = Calcular_Fadiga(Ciclo_Fadiga[i])
    
#Lê valores de Tensão Média para o caso do PSD

csv = Extrair_Dados('PSD\Media_Tensao')

D_PSD = Tamanho_Array(csv)

tensao_media = csv[:,3]

#Lê os valores de tensão para 1 Sigma do PSD

tensao_1sigma = np.zeros((2,D_PSD))

csv = Extrair_Dados('PSD\Fechado\Tensao_Equivalente_1_Sigma')

tensao_1sigma[0] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Tensao_Equivalente_1_Sigma')

tensao_1sigma[1] = csv[:,3]

#Corrigir tensao alternada para tensao media por diferentes criterios

tensao_1sigma_c = np.zeros((2,D_PSD))

#Correcao por Goodman R=0

#tensao_1sigma_c = ( tensao_1sigma * Tensao_ultima / ( ( tensao_1sigma - tensao_media ) + Tensao_ultima) ) / Fatores

#Correcao por Walker

tensao_1sigma_c = ((tensao_1sigma + tensao_media)**(1-gama) * tensao_1sigma ** gama) / Fatores

#Sem correcao

#tensao_1sigma_c = tensao_1sigma / Fatores

#Calcular frequencias

#Obtem deslocamentos
csv = Extrair_Dados('PSD\Fechado\Deslocamento_x')

deslocamento = np.zeros((6,D_PSD))

deslocamento[0] = csv[:,3]

csv = Extrair_Dados('PSD\Fechado\Deslocamento_y')

deslocamento[1] = csv[:,3]

csv = Extrair_Dados('PSD\Fechado\Deslocamento_z')

deslocamento[2] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Deslocamento_x')

deslocamento[3] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Deslocamento_y')

deslocamento[4] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Deslocamento_z')

deslocamento[5] = csv[:,3]

#Obtem velocidades

csv = Extrair_Dados('PSD\Fechado\Velocidade_x')

velocidade = np.ones((6,D_PSD))
velocidade[0] = csv[:,3]

csv = Extrair_Dados('PSD\Fechado\Velocidade_y')

velocidade[1] = csv[:,3]

csv = Extrair_Dados('PSD\Fechado\Velocidade_z')

velocidade[2] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Velocidade_x')

velocidade[3] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Velocidade_y')

velocidade[4] = csv[:,3]

csv = Extrair_Dados('PSD\Aberto\Velocidade_z')

velocidade[5] = csv[:,3]

#Calcula frequencias

frequencias = np.zeros((8,D_PSD))

for m in range(0,D_PSD,1):
    for n in range(0,6,1):
        if deslocamento[n,m] != 0:
            frequencias[n,m] = velocidade[n,m] / deslocamento[n,m]
        else:
            frequencias[n,m] = 0

for m in range(0,D_PSD,1):
    for n in range(6,8,1):
        if n == 6:
            frequencias[n,m] = np.max(frequencias[0:2,m])
        else:
            frequencias[n,m] = np.max(frequencias[3:5,m])
       
#Calcula o numero de ciclos

tempo_ensaio = 120

fator_aberto_fechado = 0.5

tempo_ensaio_aberto = tempo_ensaio * fator_aberto_fechado

tempo_ensaio_fechado = tempo_ensaio * (1 - fator_aberto_fechado)

Ciclos = np.zeros((8,D_PSD))

Ciclos[0] = frequencias[3] * tempo_ensaio_fechado
Ciclos[7] = frequencias[6] * tempo_ensaio_aberto 

#Calcula numero de ciclos em cada intensidade de vibração

Ciclos[1] = Ciclos[0]*.683
Ciclos[2] = Ciclos[0]*.271
Ciclos[3] = Ciclos[0]*.0433

Ciclos[4] = Ciclos[1]*.683
Ciclos[5] = Ciclos[1]*.271
Ciclos[6] = Ciclos[1]*.0433

#Vida para cada nivel de tensao

Ciclos_maximos = np.zeros((8,D_PSD))

Ciclos_maximos[1] = (1 * tensao_1sigma_c[0] / a) ** (1/b)
Ciclos_maximos[2] = (2 * tensao_1sigma_c[0] / a) ** (1/b)
Ciclos_maximos[3] = (3 * tensao_1sigma_c[0] / a) ** (1/b)

Ciclos_maximos[4] = (1 * tensao_1sigma_c[1] / a) ** (1/b)
Ciclos_maximos[5] = (2 * tensao_1sigma_c[1] / a) ** (1/b)
Ciclos_maximos[6] = (3 * tensao_1sigma_c[1] / a) ** (1/b)

#Calculo dos danos

Dano = np.zeros((8,D_PSD))

for i in range(1,7,1):
    Dano[i] = Ciclos[i]/Ciclos_maximos[i]
    

#Dano Total PSD

Dano[0] = Dano[1] + Dano[2] + Dano[3]
Dano[7] = Dano[4] + Dano[5] + Dano[6]

Dano_Total_PSD = np.transpose(np.array((Dano[0] + Dano[7]),ndmin=2))

for i in range(0,D_PSD,1):
    if Dano_Total_PSD[i] < (1/N_infinito):
        Dano_Total_PSD[i] = 0

Vida_PSD = np.ones((D_PSD,1))
Tensao_Equivalente_PSD = np.ones((D_PSD,1))
Coef_Seg_PSD = np.ones((D_PSD,1))


for i in range(0,D_PSD,1):
    if Dano_Total_PSD[i] == 0:
        
        Vida_PSD[i] = N_infinito
        Tensao_Equivalente_PSD[i] = 0
        Coef_Seg_PSD[i] = 5
        
    else:
        
        Vida_PSD[i] = 1/Dano_Total_PSD[i]
        Tensao_Equivalente_PSD[i] = a * Vida_PSD[i] ** b
        
        if Vida_PSD[i] / Expectativa_vida < 5:
            
            Coef_Seg_PSD[i] = Vida_PSD[i] / Expectativa_vida
            
        else:
            
            Coef_Seg_PSD[i] = 5

#Reune os resultados de PSD em um array

Resultados_Fadiga_PSD = np.concatenate((Dano_Total_PSD,Vida_PSD,Tensao_Equivalente_PSD,Coef_Seg_PSD), axis =1)
    
#Juntar os resultados de PSD e fadiga por ciclos
        
Dano_Total = np.transpose(np.array((Resultados_Fadiga_Ciclos[:,0] + Resultados_Fadiga_PSD[:,0]),ndmin=2))

Vida_Total = np.ones((D,1))
Tensao_Equivalente_Total = np.ones((D,1))
Coef_Seg_Total = np.ones((D,1))

for i in range(0,D,1):
    if Dano_Total[i] < (1 / N_infinito):
        
        Vida_Total[i] = N_infinito
        Tensao_Equivalente_Total[i] = 0
        Coef_Seg_Total[i] = 5
        
    else:
        
        Vida_Total[i] = 1/Dano_Total[i]
        Tensao_Equivalente_Total[i] = a * Vida_Total[i] ** b
        
        if Vida_Total[i] / Expectativa_vida < 5:
            
            Coef_Seg_Total[i] = Vida_Total[i] / Expectativa_vida
            
        else:
            
            Coef_Seg_Total[i] = 5
            
        

#Reune os resultados totais em um array

Resultados_Fadiga = np.concatenate((Dano_Total,Vida_Total,Tensao_Equivalente_Total,Coef_Seg_Total), axis =1)
        
#Escrever arquivo com os resultados
    
Resultados_Fadiga_coord = np.concatenate((Desmontagem[:,0:3],Resultados_Fadiga_Ciclos,Resultados_Fadiga_PSD,Resultados_Fadiga), axis =1)

np.savetxt("Dados\Fadiga\Proposta\Resultados\Resultados_Ciclos_mais_PSD_Walker.csv", Resultados_Fadiga_coord, delimiter=";",header="X;Y;Z;Dano Ciclos;Vida Ciclos;Tensao Eq Ciclos;Coef Seg Ciclos;Dano PSD;Vida PSD;Tensao Eq PSD;Coef Seg PSD;Dano Total;Vida Total;Tensao Eq Total;Coef Seg Total")