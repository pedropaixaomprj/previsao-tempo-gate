# Autor: Pedro Paixão
# Data: 19/09/2024

##################################################################################
############################## IMPORTAÇÃO DE PACOTES #############################
##################################################################################

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os

filePath=os.path.expanduser('~')+'/filesUpload/'

##################################################################################
############################### LEITURA DE ARQUIVOS ##############################
##################################################################################

admiss = pd.read_csv(filePath+'mvw_admiss.csv')
prod = pd.read_csv(filePath+'mvw_prod.csv')

admiss['ENTRADA'] = pd.to_datetime(admiss['ENTRADA'])
prod['ENTRADA'] = pd.to_datetime(prod['ENTRADA'])
prod['ADMISSAO'] = pd.to_datetime(prod['ADMISSAO'])
prod['DISTRIBUICAO'] = pd.to_datetime(prod['DISTRIBUICAO'])
prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])

##################################################################################
############################# DEFINE PARÂMETROS ##################################
##################################################################################

# Define percentil que será utilizado para obter estatística do tempo entre distribuição e criação da IT para cada tema
percentil_tempo = 0.7

# Define dias a considerar no filtro do tema
dias_a_considerar = 180

##################################################################################
####### OBTÉM ESTATÍSTICAS SOBRE TEMPO DE ADMISSÃO, DISTRIBUIÇÃO E CRIAÇÃO #######
##################################################################################

# Essa seção tem como objetivo obter estatísticas sobre o tempo de admissão, distribuição e criação
# de ITs. Consideramos os dados de ITs produzidas nos últimos 90 dias. Esses valores serão utilizados
# para prever o tempo total de produção de SATs na fase de admissibilidade, já que ainda não teremos
# informação sobre o tema nesse momento.

data_limiar = datetime.now() - timedelta(days=dias_a_considerar)

dias_admiss = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_ADMISS
dias_fila = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_FILA
dias_ate_distrib = (pd.to_datetime(prod.loc[prod['CRIACAO_IT'] > data_limiar].DISTRIBUICAO) - pd.to_datetime(prod.loc[prod['CRIACAO_IT'] > data_limiar].ENTRADA)).dt.days
dias_prod_it = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_PROD_IT
dias_total = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_TOTAL

# Gera dataframe com percentis 
percentis_possiveis = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
percentis_dias = pd.DataFrame({'Percentil': percentis_possiveis,
             'DIAS_ADMISS': dias_admiss.quantile(percentis_possiveis),
             'DIAS_FILA': dias_fila.quantile(percentis_possiveis),
             'DIAS_ATE_DISTRIB': dias_ate_distrib.quantile(percentis_possiveis),
             'DIAS_PROD': dias_prod_it.quantile(percentis_possiveis),
             'DIAS_TOTAL': dias_total.quantile(percentis_possiveis)})

##################################################################################
########################### GERA PREVISÕES DE TEMPO ##############################
##################################################################################

percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS']


# Atenção! Aqui é a previsão ATÉ a distribuição, não a previsão do tempo de fila

admiss['PREVISAO_ADMIT'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS'].values[0]
# admiss['PREVISAO_DISTRIB'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_FILA'].values[0]
admiss['PREVISAO_DISTRIB'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ATE_DISTRIB'].values[0]
admiss['PREVISAO_PROD'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_TOTAL'].values[0]

admiss['PREVISAO_ADMIT_DATA'] = pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS'].values[0])).dt.strftime('%Y-%m-%d')
admiss['PREVISAO_DISTRIB_DATA'] = pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ATE_DISTRIB'].values[0])).dt.strftime('%Y-%m-%d')
admiss['PREVISAO_PROD_DATA'] = pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_TOTAL'].values[0])).dt.strftime('%Y-%m-%d')

# Outra ideia para previsão de distribuição:
# admiss['PREVISAO_DISTRIB_DATA'] = pd.to_datetime(pd.to_datetime(admiss['PREVISAO_ADMIT_DATA']) + timedelta(days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_FILA'].values[0])).dt.strftime('%Y-%m-%d')

##################################################################################################
######################################  EXPORTA DADOS ############################################
##################################################################################################

admiss.to_csv(filePath+'previsao_tempo_admiss.csv', index=False)