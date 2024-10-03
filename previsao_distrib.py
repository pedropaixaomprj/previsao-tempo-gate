# Autor: Pedro Paixão
# Data: 19/09/2024

##################################################################################
############################## IMPORTAÇÃO DE PACOTES #############################
##################################################################################

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import math
import os

filePath=os.path.expanduser('~')+'/filesUpload/'

##################################################################################
############################### LEITURA DE ARQUIVOS ##############################
##################################################################################

# Aqui vou ter que atualizar para ler direto do banco... Ainda não consegui fazer isso pelo Spoon porque não está aceitando o código em Python.

distrib = pd.read_csv(filePath+'mvw_distrib.csv')
prod = pd.read_csv(filePath+'mvw_prod.csv')
tps = pd.read_csv(filePath+'funcs_atual.csv', sep =';')
feriados = pd.read_csv(filePath+'feriados_escalas_fds.csv')
afasts = pd.read_csv(filePath+'afastamentos_funcs.csv')

# Muda tipos para "datetime" (será necessário para realizar algumas transformações futuras)
afasts.loc[:, 'INICIO'] = pd.to_datetime(afasts.loc[:, 'INICIO'])
afasts.loc[:, 'FIM'] = pd.to_datetime(afasts.loc[:, 'FIM'])

distrib['ENTRADA'] = pd.to_datetime(distrib['ENTRADA'])
distrib['DISTRIBUICAO'] = pd.to_datetime(distrib['DISTRIBUICAO'])
prod['ENTRADA'] = pd.to_datetime(prod['ENTRADA'])
prod['ADMISSAO'] = pd.to_datetime(prod['ADMISSAO'])
prod['DISTRIBUICAO'] = pd.to_datetime(prod['DISTRIBUICAO'])
prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])

feriados['DATA'] = pd.to_datetime(feriados['DATA'], format='%d/%m/%Y')

temas = sorted(pd.Series([item for sublist in prod['TEMAS'].str.split(',') for item in sublist]).unique())
temas = [tema for tema in temas if tema not in ['PSICOLOGIA', '###   SEM TEMA   ###']]

##################################################################################
############################# DEFINE PARÂMETROS ##################################
##################################################################################

# Define percentil que será utilizado para obter estatística do tempo entre distribuição e criação da IT para cada tema
percentil_tempo = 0.70

# Define dias a considerar no filtro do tema
dias_a_considerar = 180

##################################################################################
####### OBTÉM ESTATÍSTICAS SOBRE TEMPO DE ADMISSÃO, DISTRIBUIÇÃO E CRIAÇÃO #######
##################################################################################

# Tempo de escrita se refere ao tempo entre distribuição e produção
tempo_escrita_temas_nACP = {}
tempo_escrita_temas_ACP = {}

for tema in temas:
  
    prod_tema_recente = prod.loc[(prod.TEMAS.str.contains(tema)) & (pd.to_datetime(prod.CRIACAO_IT) >= datetime.today() - timedelta(days=dias_a_considerar))]
 
    if len(prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 5:  # Ao menos 5 não-ACPs produzidas nos últimos 180 dias

        # Cálculo usando apenas produções recentes do tema
        tempo_escrita_temas_nACP[tema] = prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)
        
    else:
        # Cálculo histórico, considerando o número de funcionários atuais no tema
        tempo_escrita_temas_nACP[tema] = prod.loc[prod.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

    # Geralmente tem poucas ACPs, então usamos o tempo histórico de produção média
    # Cálculo histórico, considerando o número de funcionários atuais no tema

    if len(prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO == 'PRAZO PROCESSUAL']) > 5:  # Ao menos 5 ACPs produzidas nos últimos 180 dias
        
        # Cálculo usando apenas produções recentes do tema
        tempo_escrita_temas_ACP[tema] = prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)
        
    else:

        # Se não tem ACP suficiente, usa os dados das SATs gerais para o tema
        
        prod_tema = prod.loc[prod.TEMAS.str.contains(tema)]

        # Se tem pelo menos 5 ACPs historicamente no tema, usa esses dados
        
        if len(prod_tema.loc[prod_tema.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 5:
            tempo_escrita_temas_ACP[tema] = prod_tema.loc[prod_tema.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

        # Se não, usa a média de ACPs no geral
        
        else:
            tempo_escrita_temas_ACP[tema] = prod.loc[prod.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

##################################################################################
########################### GERA PREVISÕES DE TEMPO ##############################
##################################################################################

# Gera coluna com tempo de previsão (em dias)

distrib['PREVISAO_PROD'] = np.where(distrib['TIPO_PRAZO'] != 'PRAZO PROCESSUAL',
                                    np.ceil(distrib['TEMA'].map(tempo_escrita_temas_nACP)),  # Map TEMA values to the dictionary for non-ACP
                                    np.ceil(distrib['TEMA'].map(tempo_escrita_temas_ACP)))   # Map TEMA values to the dictionary for ACP

distrib['PREVISAO_PROD'] = distrib['PREVISAO_PROD'].fillna(np.nan)

# Gera coluna com data final de produção aplicando função a todas as linhas

def calcula_data_previsao(row):
    if pd.isna(row['TEMA']):
        # Se o tema for nan, retorna nan
        return pd.NaT  # Return a missing value if 'TEMA' is NaN
    if row['TIPO_PRAZO'] != 'PRAZO PROCESSUAL':
        # Se não for ACP, retorna o tempo esperado para não-ACPs no tema
        return pd.to_datetime(row['DISTRIBUICAO']) + timedelta(days=tempo_escrita_temas_nACP[row['TEMA']])
    else:
        # Se for ACP, retorna o tempo esperado para ACPs no tema
        return pd.to_datetime(row['DISTRIBUICAO']) + timedelta(days=tempo_escrita_temas_ACP[row['TEMA']])

distrib['PREVISAO_PROD_DATA'] = distrib.apply(calcula_data_previsao, axis=1)

# Função para transformar finais de semana em segundas
def proximo_dia_util(date):
    if pd.isna(date):  # Se é nan
        return date
    if date.weekday() >= 5:  # Se é sábado (5) ou domingo (6)
        return date + pd.offsets.BDay(1)
    return date  # Se já é um dia útil

# Aplica a função e transforma em string
distrib['PREVISAO_PROD_DATA'] = distrib['PREVISAO_PROD_DATA'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d')

##################################################################################################
######################################  EXPORTA DADOS ############################################
##################################################################################################

distrib.to_csv(filePath+'previsao_tempo_distrib.csv', index=False)