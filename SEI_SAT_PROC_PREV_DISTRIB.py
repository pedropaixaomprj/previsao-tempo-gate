import psycopg2
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

conn = psycopg2.connect(
    host="h-pgsql01.pgj.rj.gov.br",
    database="gate",
    user="gate",
    password="gatehom2020",
    port="5432"  # Default PostgreSQL port
)

conn.autocommit = True

cursor = conn.cursor()

query_prod = '''
SELECT "ENTRADA"
     , "ADMISSAO"
     , "DISTRIBUICAO"
     , "CRIACAO_IT"
     , "TIPO_PRAZO"
     , "DIAS_ADMISS"
     , "DIAS_FILA"
     , "DIAS_PROD_IT"
     , "DIAS_TOTAL"
     , "TEMAS"
FROM stage."MVW_SEI_SAT_GESTAO_ACERVO_PROD" prod
'''

query_distrib = '''
select "SEI",
		"ID_PROTOCOLO",
		"DISTRIBUICAO",
		"TIPO_PRAZO",
		"TEMA"
from stage."MVW_SEI_SAT_GESTAO_ACERVO_DISTRIB"'''

prod = pd.read_sql(query_prod, conn)

distrib = pd.read_sql(query_distrib, conn)

prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])

##################################################################################
############################# DEFINE PARÂMETROS ##################################
##################################################################################

# Lista temas a serem incluídos
temas = sorted(pd.Series([item for sublist in prod['TEMAS'].str.split(',') for item in sublist]).unique())
temas = [tema for tema in temas if tema not in ['PSICOLOGIA', '###   SEM TEMA   ###']]

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

    prod_tema_recente = prod.loc[(prod.TEMAS.str.contains(tema)) & (
                pd.to_datetime(prod.CRIACAO_IT) >= datetime.today() - timedelta(days=dias_a_considerar))]

    if len(prod_tema_recente.loc[
               prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 5:  # Ao menos 5 não-ACPs produzidas nos últimos 180 dias

        # Cálculo usando apenas produções recentes do tema
        tempo_escrita_temas_nACP[tema] = prod_tema_recente.loc[
            prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

    else:
        # Cálculo histórico, considerando o número de funcionários atuais no tema
        tempo_escrita_temas_nACP[tema] = prod.loc[prod.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(
            percentil_tempo)

    # Geralmente tem poucas ACPs, então usamos o tempo histórico de produção média
    # Cálculo histórico, considerando o número de funcionários atuais no tema

    if len(prod_tema_recente.loc[
               prod_tema_recente.TIPO_PRAZO == 'PRAZO PROCESSUAL']) > 5:  # Ao menos 5 ACPs produzidas nos últimos 180 dias

        # Cálculo usando apenas produções recentes do tema
        tempo_escrita_temas_ACP[tema] = prod_tema_recente.loc[
            prod_tema_recente.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

    else:

        # Se não tem ACP suficiente, usa os dados das SATs gerais para o tema

        prod_tema = prod.loc[prod.TEMAS.str.contains(tema)]

        # Se tem pelo menos 5 ACPs historicamente no tema, usa esses dados

        if len(prod_tema.loc[prod_tema.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 5:
            tempo_escrita_temas_ACP[tema] = prod_tema.loc[
                prod_tema.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

        # Se não, usa a média de ACPs no geral

        else:
            tempo_escrita_temas_ACP[tema] = prod.loc[prod.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(
                percentil_tempo)

##################################################################################
########################### GERA PREVISÕES DE TEMPO ##############################
##################################################################################

# Gera coluna com tempo de previsão (em dias)

distrib['PREVISAO_PROD'] = np.where(distrib['TIPO_PRAZO'] != 'PRAZO PROCESSUAL',
                                    np.ceil(distrib['TEMA'].map(tempo_escrita_temas_nACP)),
                                    # Map TEMA values to the dictionary for non-ACP
                                    np.ceil(distrib['TEMA'].map(
                                        tempo_escrita_temas_ACP)))  # Map TEMA values to the dictionary for ACP

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
distrib['PREVISAO_PROD_DATA'] = pd.to_datetime(
    distrib['PREVISAO_PROD_DATA'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d'))

###### LIMPA stage."SEI_SAT_PROC_PREV_DISTRIB" e insere as previsões



truncate_query = 'truncate table stage."SEI_SAT_PROC_PREV_DISTRIB"'

# Execute the command
cursor.execute(truncate_query)

distrib['PREVISAO_PROD'] = distrib['PREVISAO_PROD'].replace({np.nan: 10000})
distrib['PREVISAO_PROD_DATA'] = distrib['PREVISAO_PROD_DATA'].replace({np.nan: '9999-01-01'})

# Insert DataFrame rows one by one
for i,row in distrib.iterrows():
    qins = '''
INSERT INTO stage."SEI_SAT_PROC_PREV_DISTRIB"
("ID_PROTOCOLO", "PREVISAO_PROD", "PREVISAO_PROD_DATA")
VALUES(%i, %i, '%s');
        ''' % (row['ID_PROTOCOLO']
           , row['PREVISAO_PROD']
           , row['PREVISAO_PROD_DATA']
        )
    cursor.execute(qins)

# Close the cursor and connection
cursor.close()
conn.close()
