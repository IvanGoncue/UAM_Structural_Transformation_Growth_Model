# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:32:21 2023

@author: GonCue
"""
#---------------------------------FUNCIONES MODELO------------------------------------
def k0(x):#la0 dada
    ka=alfN*(1-alfA)*x/(alfA*(1-alfN))
    MPKn=(1-alfN)*Bn*(An/ (Bn*x) )** (alfN)
    MPKa = (1-alfA)*Ba*(Aa/(Ba*ka))**(alfA)
    pa=MPKn/MPKa
    ya=la0*Aa**(alfA)*(Ba*ka)**(1-alfA)
    c=(ya-c_a)*(pa)**(1-omgA)/omgA
    yn=(1-la0)*An**(alfN)*(Bn*x)**(1-alfN)
    fin=(1-s)*yn-(1-omgA)*(pa)**(omgA)*c
    return fin

def kk(x,valor):
    kk=valor
    ka=alfN*(1-alfA)*x/(alfA*(1-alfN))
    MPKn=(1-alfN)*Bn*(An/ (Bn*x) )** (alfN)
    MPKa=(1-alfA)*Ba*(Aa/(Ba*ka))**(alfA)
    pa=MPKn/MPKa
    yaa=Aa**(alfA)*(Ba*ka)**(1-alfA)
    ybb=An**(alfN)*(Bn*x)**(1-alfN)
    num=ybb*(1-s)+c_a*(1-omgA)*pa/omgA
    den=ybb*(1-s)+yaa*(1-omgA)*pa/omgA
    la=num/den
    fin=la*ka+(1-la)*x-kk
    return fin
#----------------------------------OPTIMIZACIÓN--------------------------------------------
from scipy.optimize import fsolve
alfA,alfN,omgA,c_a=0.47,0.687,0.005,0.35
s,la0=0.15,0.8
Aa,An,Ba,Bn=1,1,0.5,0.5
gAa,gAn,gBa,gBn,delta,start=0.010,0.012,0,0,0.05,0.5


k,kn,la,ka,ya,yn,t,pa=([]for i in range (8))
kn0=fsolve(k0,start)
ka0=alfN*(1-alfA)*kn0/(alfA*(1-alfN))
k0=la0*ka0+(1-la0)*kn0
ya0=la0*Aa**(alfA)*(Ba*ka0)**(1-alfA)
yn0=(1-la0)*An**(alfN)*(Bn*kn0)**(1-alfN)
k1=s*yn0+(1-delta)*k0
MPKn0=(1-alfN)*Bn*(An/ (Bn*kn0) )** (alfN)
MPKa0=(1-alfA)*Ba*(Aa/(Ba*ka0))**(alfA)
pa0=MPKn0/MPKa0

k.append(k0),k.append(k1)
kn.append(kn0),ka.append(ka0)
la.append(la0),pa.append(pa0)
ya.append(ya0),yn.append(yn0)
t.append(0)

#-----------------------------------"N" PERIODOS--------------------------------------------
years=range(500)
for i in years:
    t.append(i+1)
    Ba=Ba*(1+gBa)
    Bn=Bn*(1+gBn)
    Aa=Aa*(1+gAa)
    An=An*(1+gAn)
    kn1=fsolve(kk,kn0,k1)
    ka1=alfN*(1-alfA)*kn1/(alfA*(1-alfN))
    MPKn1=(1-alfN)*Bn*(An/(Bn*kn1) )**(alfN)
    MPKa1=(1-alfA)*Ba*(Aa/(Ba*ka1))**(alfA)
    pa1=MPKn1/MPKa1
    yaa1=Aa**(alfA)*(Ba*ka1)**(1-alfA)
    ybb1=An**(alfN)*(Bn*kn1)**(1-alfN)
    num=ybb1*(1-s)+c_a*(1-omgA)*pa1/omgA
    den=ybb1*(1-s)+yaa1*(1-omgA)*pa1/omgA
    la1=num/den
    ya1=yaa1*la1
    yn1=ybb1*(1-la1)
    k2=s*yn1+(1-delta)*k1
    k1=k2
    k.append(k1)
    kn.append(kn1),ka.append(ka1)
    la.append(la1),pa.append(pa1)
    ya.append(ya1),yn.append(yn1)


#----------------------------------GRÁFICAMENTE--------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
#shares en output
myn=np.array(yn)
mya=np.array(ya)
mpa=np.array(pa)
sharagr=mpa*mya/(mpa*mya+myn)
sharnagr=myn/(mpa*mya+myn)

#GRÁFICO 1---- OUTPUT SHARES
# Crea la figura con tres subplots
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20,5))
fig.set_dpi(300)
# Primer subplot: PESOS SECTORIALES RELATIVOS/t
axs[0].plot(sharagr, label="Share Agric")
axs[0].plot(sharnagr, label="Share No agric")
axs[0].legend(loc="center right")
axs[0].set_xlabel("Tiempo")
axs[0].set_ylabel("Valores")
axs[0].set_title("PESOS SECTORIALES RELATIVOS/t")

# Segundo subplot: PESOS SECTORIALES RELATIVOS/yn
axs[1].plot(sharagr, label="pa*ya/(pa*ya+yn)")
axs[1].plot(sharnagr, label="yn/(pa*ya+yn)")
axs[1].legend(loc="center right")
axs[1].set_xlabel("Valores yn(0), yn(60) e yn(140)")
axs[1].set_ylabel("Valores")
axs[1].set_title("PESOS SECTORIALES RELATIVOS/yn")
t=[0, 60, 140]
axs[1].set_xticks(t)
axs[1].set_xticklabels([yn[i] for i in t])

# Tercer subplot: PESOS SECTORIALES RELATIVOS/ya
axs[2].plot(sharagr, label="pa*ya/(pa*ya+yn)")
axs[2].plot(sharnagr, label="yn/(pa*ya+yn)")
axs[2].legend(loc="center right")
axs[2].set_xlabel("Valores ya(0), ya(60) e ya(140)")
axs[2].set_ylabel("Valores")
axs[2].set_title("PESOS SECTORIALES RELATIVOS/ya")
t=[0, 60, 140]
axs[2].set_xticks(t)
axs[2].set_xticklabels([ya[i] for i in t])
#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
# Muestra la figura completa
plt.show()

#GRÁFICO 2----Todas las variables
fig, axa = plt.subplots()
fig.set_dpi(300)
axa.plot(k, label="k")
axa.plot(kn, label="kn")
axa.plot(ka, label="ka")
axa.plot(la, label="la")
axa.plot((1-np.array(la)), label="ln")
axa.plot(pa, label="pa")
axa.plot(ya, label="ya")
axa.plot(yn, label="yn")
axa.legend(loc="upper left")
axa.set_xlabel("Tiempo")
axa.set_ylabel("Valores")
axa.set_title("EVOLUCIÓN VARIABLES")
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
plt.show()

#GRÁFICO 3---k
fig, axa = plt.subplots()
fig.set_dpi(300)
axa.plot(k, label="k")
axa.legend(loc="upper left")
axa.set_xlabel("Tiempo")
axa.set_ylabel("Valores")
axa.set_title("EVOLUCIÓN k")
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
plt.show()

#GRÁFICO 4----EVOLUCIÓN CAPITAL POR SECTORES (kn/ka)
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20,5),sharey=True)
fig.set_dpi(300)#shrey=true pone la misma escala eje y

# Primer subplot: CAPITAL NO AGRÍCOLA
axs[0].plot(kn, label="kn")
axs[0].legend(loc="center right")
axs[0].set_xlabel("Tiempo")
axs[0].set_ylabel("Valores")
axs[0].set_title("EVOLUCIÓN CAPITAL NO AGRÍCOLA")

# Segundo subplot: CAPITAL AGRÍCOLA
axs[1].plot(ka, label="ka")
axs[1].legend(loc="center right")
axs[1].set_xlabel("Tiempo")
axs[1].set_ylabel("Valores")
axs[1].set_title("EVOLUCIÓN CAPITAL AGRÍCOLA")


#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
# Muestra la figura completa
plt.show()

#GRÁFICO 5----EVOLUCIÓN PROPORCIÓN DE OCUPACIÓN POR SECTORES (ln/la)
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20,5))
fig.set_dpi(300)

# Primer subplot: PROPORCIÓN OCUPACIÓN NO AGRÍCOLA
axs[0].plot((1-np.array(la)), label="ln=(1-la)")
axs[0].legend(loc="center right")
axs[0].set_xlabel("Tiempo")
axs[0].set_ylabel("Valores")
axs[0].set_title("EVOLUCIÓN PROPORCIÓN DE OCUPACIÓN NO AGRÍCOLA")

# Segundo subplot: PROPORCIÓN OCUPACIÓN AGRÍCOLA
axs[1].plot(la, label="la")
axs[1].legend(loc="center right")
axs[1].set_xlabel("Tiempo")
axs[1].set_ylabel("Valores")
axs[1].set_title("EVOLUCIÓN PROPORCIÓN DE OCUPACIÓN AGRÍCOLA")


#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
# Muestra la figura completa
plt.show()

#GRÁFICO 6---- EVOLUCIÓN PRODUCCIÓN POR SECTORES (yn/ya)
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20,5))
fig.set_dpi(300)

# Primer subplot: PRODUCCIÓN NO AGRÍCOLA
axs[0].plot(yn, label="yn")
axs[0].legend(loc="center right")
axs[0].set_xlabel("Tiempo")
axs[0].set_ylabel("Valores")
axs[0].set_title("EVOLUCIÓN PRODUCCIÓN NO AGRÍCOLA")

# Segundo subplot: PRODUCCIÓN AGRÍCOLA
axs[1].plot(ya, label="ya")
axs[1].legend(loc="center right")
axs[1].set_xlabel("Tiempo")
axs[1].set_ylabel("Valores")
axs[1].set_title("EVOLUCIÓN PRODUCCIÓN AGRÍCOLA")
axs[1].set_ylim(axs[0].get_ylim()) # establecer los mismos límites de eje

#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
# Muestra la figura completa
plt.show()

#GRÁFICO 7---PRECIO AGRÍCOLA (pa)
fig, axa = plt.subplots()
fig.set_dpi(300)
axa.plot(pa, label="pa")
axa.legend(loc="upper left")
axa.set_xlabel("Tiempo")
axa.set_ylabel("Valores")
axa.set_title("EVOLUCIÓN pa")
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')
plt.show()


#-------------------------------------ECONOMÍA USA-----------------------------------

import pandas as pd
datos=pd.read_csv('datos.csv', header=0)#header=0 significa que la primera fila es el encabezado
datos = datos.iloc[25:].reset_index(drop=True)
#para deshacerme de los datos previos a 1925
print(datos)


la_real=datos["la"]
vaava=datos["vaava"]
kakn=datos["kakn"]
papn=datos["papn"]
YiYUS1925=datos["YiYUS1925"]
year=datos['year']
    
#extraer los datos de una de las columnas
#print(datos['countryname'])#se pone el nombre de la columna entrecomillado y en corchetes
#print(datos.iloc[0:3])#este programa me da las 3 primeras filas
#antes era .ix, desde pandas 0.20 es loc o iloc
#print(datos.sort_values(by='year', ascending=False))#ordeno los datos, en este caso por año de < a >
#si quitamos ascending=False los ordena de menor a mayor

# Crear una figura y un conjunto de subplots
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10, 7), gridspec_kw={'hspace': 0.8, 'wspace': 0.2})
fig.suptitle('COMPARACIÓN DEL MODELO CON LA ECONOMÍA EE.UU 1925-2010', fontsize=16)
plt.tight_layout()#para que los gráficos no se sobrepongan

#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')

# Dibujar cada subplot con los datos correspondientes
fig.set_dpi(300)

#IMPORTANTE!!!!!!!!¡¡¡¡¡¡¡
#posición l_a=0.25-->período x del modelo, como l_a en EEUU 1925 es 0.25
#la simulación la cogemos desde el periodo en el que la l_a es la misma que en US. 

#posición la=0,25

valor_buscado = 0.246994 #l_a US 1925
tolerancia = 0.03

posicion = None
for i, valor in enumerate(la):
    if np.isclose(valor, valor_buscado, rtol=tolerancia):
        posicion = i
        break

if posicion is not None:
    print(f"El valor {valor_buscado} se encontró en la lista en la posición {posicion}.")
else:
    print(f"No se encontró el valor {valor_buscado} en la lista.")


#PRIMER SUBPLOT
axs[0, 0].plot(la_real, label='Datos Reales-la', alpha=0.8)
axs[0, 0].plot(la[posicion:141], color='red', label='Datos Modelo-la',alpha=0.5)
axs[0, 0].set_title('PROPORCIÓN DE OCUPACIÓN AGRÍCOLA')
t=[0,20,40,60,80]
axs[0, 0].set_xticks(t)
axs[0, 0].set_xticklabels([year[i] for i in t])


#SEGUNDO SUBPLOT
ln=1-np.array(la)
axs[0, 1].plot(1-la_real, label='Datos Reales-ln', alpha=0.8)
axs[0, 1].plot(ln[posicion:141], label='Datos Modelo-ln', color='red',alpha=0.5)
axs[0, 1].set_title('PROPORCIÓN DE OCUPACIÓN NO AGRÍCOLA')
t=[0,20,40,60,80]
axs[0, 1].set_xticks(t)
axs[0, 1].set_xticklabels([year[i] for i in t])

#TERCER SUBPLOT
kaknmodel=np.array(ka)/np.array(kn)
axs[1, 0].plot(kakn, label='Datos Reales-ka/kn',alpha=0.8)
axs[1, 0].plot(kaknmodel[posicion:141], label='Datos Modelo-ka/kn',color='red',alpha=0.5)
axs[1, 0].set_title('RATIO DE LOS CAPITALES SECTORIALES')
t=[0,20,40,60,80]
axs[1, 0].set_xticks(t)
axs[1, 0].set_xticklabels([year[i] for i in t])

#CUARTO SUBPLOT
pamodel=np.array(pa)
axs[1, 1].plot(papn, label='Datos Reales-pa/pn',alpha=0.8)
#axs[1, 1].plot((pamodel[:88]/pamodel[0]), label='Datos Modelo-pa/pn', color='red',alpha=0.5 )#dividir p[0] para 1er valor=1
axs[1, 1].plot((pamodel[posicion:141]/pamodel[posicion]), label='Datos Modelo-pa/pn', color='red',alpha=0.5 )
axs[1, 1].set_title('RATIO DE LOS PRECIOS SECTORIALES')
t=[0,20,40,60,80]
axs[1, 1].set_xticks(t)
axs[1, 1].set_xticklabels([year[i] for i in t])


#QUINTO SUBPLOT
gdp=np.array(pa)*np.array(ya)+np.array(yn)
gdpindex=gdp/gdp[posicion]
axs[2, 0].plot(YiYUS1925, label='Datos Reales-y',alpha=0.8)
axs[2, 0].plot(gdpindex[posicion:142], label='Datos Modelo-y',color='red',alpha=0.5 )
axs[2, 0].set_title('PIB per cápita index la=0,25')
t=[0,20,40,60,80]
axs[2, 0].set_xticks(t)
axs[2, 0].set_xticklabels([year[i] for i in t])
# evolución de la producción en el índice de posición de la proporción 
#de ocupación agrícola, que es igual a 0,25 en los datos reales

   
# SEXTO SUBPLOT
axs[2, 1].plot(vaava/vaava[0], label='Datos Reales-vaagr',alpha=0.8)
axs[2, 1].plot(sharagr[posicion:141]/sharagr[posicion], label='Datos Modelo-vaagr', color='red',alpha=0.5 )
axs[2, 1].set_title('VALOR AÑADIDO SECTOR AGRÍCOLA/TOTAL')
t=[0,20,40,60,80]
axs[2, 1].set_xticks(t)
axs[2, 1].set_xticklabels([year[i] for i in t])

#-----------------------------------CUESTIÓN PA--------------------------------
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), gridspec_kw={'hspace': 0.8, 'wspace': 0.2})
fig.suptitle('RATIO DE LOS PRECIOS SECTORIALES', fontsize=16)
pamodel=np.array(pa)
axs.plot(papn,label='Precios reales since la=0.25')
axs.plot((pamodel[:88]/pamodel[0]),label='Precios modelo since la=0.8')#dividir p[0] para 1er valor=1
t=[0,20,40,60,80]
axs.set_xticks(t)
axs.set_xticklabels([year[i] for i in t])
axs.legend()
#------------------------------------CROSS-COUNTRY-----------------------------------
dat=pd.read_csv('crosscountry.csv', header=0)#header=0 significa que la primera fila es el encabezado

vaava2 = dat['vaava']
la_asger = dat['la_asger']
ln_asger = dat['ln_asger']
kakn_tot_asger = dat['kakn_tot_asger']
papn_narrow1 = dat['papn_narrow1']
yiyUS = dat['yiyUS']

print(dat)

# Crear una figura y un conjunto de subplots
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 7), gridspec_kw={'hspace': 0.5, 'wspace': 0.2})
fig.suptitle('COMPROBACIÓN EMPÍRICA DEL MODELO', fontsize=16)
plt.tight_layout()#para que los gráficos no se sobrepongan

#fondo y borrde del gráfico
fig.set_facecolor('#D3D3D3')
fig.set_edgecolor('black')

# Dibujar cada subplot con los datos correspondientes
fig.set_dpi(300)


#PRIMER SUBPLOT
axs[0, 0].plot(ln,sharagr,color='red',alpha=0.7)
axs[0, 0].plot(ln_asger,vaava2,'.',markersize=12, alpha=0.5)
axs[0, 0].set_title('VALOR AÑADIDO SECTOR AGRÍCOLA/TOTAL')
axs[0, 0].set_xlabel('Proporción de ocupación no agrícola')
# Calcular los coeficientes de la recta de regresión
m, b = np.polyfit(ln_asger, vaava2, 1)
# Dibujar la recta de regresión sobre los datos
axs[0, 0].plot(ln_asger, m * ln_asger + b, '--',color='black', alpha=0.5, linewidth=2)



#SEGUNDO SUBPLOT 
valor_buscado = 0.015 #l_a US 2010
tolerancia = 0.05

pos = None
for i, valor in enumerate(la):
    if np.isclose(valor, valor_buscado, rtol=tolerancia):
        pos = i
        break

if pos is not None:
    print(f"El valor {valor_buscado} se encontró en la lista en la posición {posicion}.")
else:
    print(f"No se encontró el valor {valor_buscado} en la lista.")

gdpp=gdp/gdp[pos]
axs[0, 1].plot(ln_asger,yiyUS,'.',alpha=0.5,markersize=12)
axs[0, 1].plot(ln[:300],gdpp[:300], color='red', alpha=0.7)
axs[0, 1].set_title('PIB per cápita')
axs[0, 1].set_xlabel('Proporción de ocupación no agrícola')
# Calcular los coeficientes de la recta de regresión
m, b = np.polyfit(ln_asger, yiyUS, 1)
# Dibujar la recta de regresión sobre los datos
#axs[0, 1].plot(ln_asger, m * ln_asger + b, '--',color='black', alpha=0.5, linewidth=2)



#TERCER SUBPLOT
axs[1, 0].plot(ln_asger,kakn_tot_asger,'.',alpha=0.5,markersize=12)
axs[1, 0].plot(ln,kaknmodel, color='red', alpha=0.7)
axs[1, 0].set_title('RATIO DE LOS CAPITALES SECTORIALES')
axs[1, 0].set_xlabel('Proporción de ocupación no agrícola')



#CUARTO SUBPLOT
axs[1, 1].plot(ln_asger,papn_narrow1,'.',alpha=0.5,markersize=12)
axs[1, 1].plot(ln,pamodel, color='red', alpha=0.7)
axs[1, 1].set_title('RATIO DE LOS PRECIOS SECTORIALES')
axs[1, 1].set_xlabel('Proporción de ocupación no agrícola')




