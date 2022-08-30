# -*- coding: utf-8 -*-
"""Extraccion Reddit.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1_jZz0T0VNtCPG62wCret539IutwZ2fPD

## Descarga de información de Reddit, mediante requests.

Se cargan los modulos necesarios para la tarea
"""

import pandas as pd # Uso de dataframes
import requests # Método para acceder a Pushshift por mediante url 
import json # Manipulación de información JSON
import csv # Para convertir las tablas finales en archivos csv y guardarla en la máquina local
import time # Convertir el tiempo UTC en tiempo GMT
import datetime # Manipulación del tiempo

# Adaptado de https://gist.github.com/dylankilkenny/3dbf6123527260165f8c5c3bc3ee331b
# Se construye la URL de Pushshift, accediendo a una pagina web donde está almacenada la información en JSON a manera de lista
def getPushshiftData(query, after, before, sub):
    # Construye la URL de Pushshift con los parámetros antes mencionados
    url = 'https://api.pushshift.io/reddit/search/submission/?title='+str(query)+'&size=1000&after='+str(after)+'&before='+str(before)+'&subreddit='+str(sub)
    # Muestra al usuario la URL 
    print(url)
    # Solicitud de la URL
    r = requests.get(url, headers = {'User-agent': 'bot_1'})
    r.raise_for_status()
    # Carga la pagina con la información JSON en una variable
    if r.status_code != 204:
      data = json.loads(r.text, strict=False,encoding='utf-8')
    #regresa la información de la variable donde está contenida
    return data['data']

# Extracción de la información específica de los datos JSON generados  
# Seleccionar los parámetros deseados de pushshift: https://pushshift.io/api-parameters/
def collectSubData(subm):
    #subData fue creado al inicio donde está toda la información para ser agregados a la variable global subStats.
    subData = list() # lista para almacenar la información
    title = subm['title'] # Titulo 
    url = subm['url'] # URL
    subreddit = subm['subreddit'] # Grupo Subreddit
    # selftext, o cuerpo del post, no siempre está presente en submissions, para evitar errores se usa try/except
    try:
      body = subm['selftext']
    except KeyError:
      body = ''
    #flairs no siempre está presente en submissions, para evitar errores se usa try/except
    try:
        flair = subm['link_flair_text']
    except KeyError:
        flair = "NaN"    
    author = subm['author'] # Autor del post
    sub_id = subm['id'] # ID del post
    score = subm['score'] # Puntuación del post
    created = datetime.datetime.fromtimestamp(subm['created_utc']) # Fecha de creación en UTC, por lo que se hace la conversión
    numComms = subm['num_comments'] # Número de comentarios
    permalink = subm['permalink'] # Linl permanente

    #Junta toda la información en una tubla y se agrega a subData
    subData.append((sub_id,title,url,subreddit,body,author,score,created,numComms,permalink,flair))
    # Crea un diccionario con la entrada de información del submission actual y almacena toda la información relacionada con el submission
    subStats[sub_id] = subData

"""## Parámetros de búsqueda"""

# Para tener el formato correcto de tiempo de la URL se puede usar la siguiente página
#https://www.unixtimestamp.com/index.php > Esto para crear el timestamp deseado
after = "1577836800" #Submissions despues del timestamp 01 Enero 2020 0:00:0 1577836800
before = "1640995200" #Submissions antes del timestamp 01 Enero 2022 0:00:0 1640995200
query = "" # Palabra clave para buscar en submissions
sub = "wallstreetbets" #Subreddit de la información

#subCount cuenta el no. del total de envíos que se recopila
subCount = 0
#subStats es el diccionario donde almacenaremos nuestros datos.
subStats = {}

data = getPushshiftData(query, after, before, sub)
# Se ejecutará hasta que se hayan recopilado todas las publicaciones, es decir, cuando la longitud de la variable de datos = 0
# desde la fecha 'posterior' hasta la fecha 'anterior'
while len(data) > 0: #La longitud de los datos es el número de submission(data[0],data[1], etc.), una vez que llega a cero (después y antes de que vars sean iguales) finaliza
    for submission in data:
        collectSubData(submission)
        subCount+=1
    # Llama getPushshiftData() Con la fecha creada del ultimo submission
    print(len(data))
    print(str(datetime.datetime.fromtimestamp(data[-1]['created_utc'])))
    #actualiza la variable 'after' a la última fecha creada del submission
    after = data[-1]['created_utc']
    # datos han cambiado debido a la nueva variable posterior proporcionada por el código anterior
    data = getPushshiftData(query, after, before, sub)
    
print(len(data))

print(str(len(subStats)) + " submissions agregados a la lista")
print("1st entrada es:")
print(list(subStats.values())[0][0][1] + "; creado por: " + str(list(subStats.values())[0][0][5]))
print("Ultima entrada es:")
print(list(subStats.values())[-1][0][1] + "; creado por: " + str(list(subStats.values())[-1][0][5]))

"""## Almacenamiento en csv"""

def updateSubs_file(location):
    upload_count = 0
    print("Ingresa el nombre del archivo agregando .csv al final ")
    #location = "/content/drive/MyDrive/Pushshift/Parts/"
    filename = input() #Nombre que puso al archivo csv
    file = location + filename
    with open(file, 'w', newline='', encoding='utf-8') as file: 
        a = csv.writer(file, delimiter=',')
        headers = ["Post ID","Title","Url","Subreddit","Body","Author","Score","Publish Date","Total No. of Comments","Permalink","Flair"]
        a.writerow(headers)
        for sub in subStats:
            a.writerow(subStats[sub][0])
            upload_count+=1
            
        print(str(upload_count) + " submissions han sido guardados en un archivo csv")

updateSubs_file(location)