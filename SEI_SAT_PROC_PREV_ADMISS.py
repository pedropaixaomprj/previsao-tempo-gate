import psycopg2
import pandas as pd
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

query_admiss = '''
select "SEI",
	"ID_PROTOCOLO",
	"ENTRADA"
from stage."MVW_SEI_SAT_GESTAO_ACERVO_ADMISS"
where "ENTRADA" is not null
'''

prod = pd.read_sql(query_prod, conn)

admiss = pd.read_sql(query_admiss, conn)

prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])

percentil_tempo = 0.7

dias_a_considerar = 180

data_limiar = datetime.now() - timedelta(days=dias_a_considerar)

dias_admiss = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_ADMISS
dias_fila = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_FILA
dias_ate_distrib = (pd.to_datetime(prod.loc[prod['CRIACAO_IT'] > data_limiar].DISTRIBUICAO) - pd.to_datetime(
    prod.loc[prod['CRIACAO_IT'] > data_limiar].ENTRADA)).dt.days
dias_prod_it = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_PROD_IT
dias_total = prod.loc[prod['CRIACAO_IT'] > data_limiar].DIAS_TOTAL

percentis_possiveis = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
percentis_dias = pd.DataFrame({'Percentil': percentis_possiveis,
                               'DIAS_ADMISS': dias_admiss.quantile(percentis_possiveis),
                               'DIAS_FILA': dias_fila.quantile(percentis_possiveis),
                               'DIAS_ATE_DISTRIB': dias_ate_distrib.quantile(percentis_possiveis),
                               'DIAS_PROD': dias_prod_it.quantile(percentis_possiveis),
                               'DIAS_TOTAL': dias_total.quantile(percentis_possiveis)})

admiss['ENTRADA'] = pd.to_datetime(admiss['ENTRADA'])
prod['ENTRADA'] = pd.to_datetime(prod['ENTRADA'])
prod['ADMISSAO'] = pd.to_datetime(prod['ADMISSAO'])
prod['DISTRIBUICAO'] = pd.to_datetime(prod['DISTRIBUICAO'])
prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])

##################################################################################
########################### GERA PREVISÕES DE TEMPO ##############################
##################################################################################

percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS']

# Atenção! Aqui é a previsão ATÉ a distribuição, não a previsão do tempo de fila

admiss['PREVISAO_ADMIT'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS'].values[0]
# admiss['PREVISAO_DISTRIB'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_FILA'].values[0]
admiss['PREVISAO_DISTRIB'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ATE_DISTRIB'].values[
    0]
admiss['PREVISAO_PROD'] = percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_TOTAL'].values[0]

admiss['PREVISAO_ADMIT_DATA'] = pd.to_datetime(pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(
    days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ADMISS'].values[0])).dt.strftime(
    '%Y-%m-%d'))
admiss['PREVISAO_DISTRIB_DATA'] = pd.to_datetime(pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(
    days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_ATE_DISTRIB'].values[0])).dt.strftime(
    '%Y-%m-%d'))
admiss['PREVISAO_PROD_DATA'] = pd.to_datetime(pd.to_datetime(pd.to_datetime(admiss['ENTRADA']) + timedelta(
    days=percentis_dias.loc[percentis_dias.Percentil == percentil_tempo, 'DIAS_TOTAL'].values[0])).dt.strftime(
    '%Y-%m-%d'))

###### LIMPA stage."SEI_SAT_PROC_PREV_ADMISS" e insere as previsões

truncate_query = 'truncate table stage."SEI_SAT_PROC_PREV_ADMISS"'

# Execute the command
cursor.execute(truncate_query)

# Insert DataFrame rows one by one
for i,row in admiss.iterrows():
    qins = '''INSERT INTO stage."SEI_SAT_PROC_PREV_ADMISS"
            ("ID_PROTOCOLO", "PREVISAO_ADMIT", "PREVISAO_DISTRIB", "PREVISAO_PROD", "PREVISAO_ADMIT_DATA",
             "PREVISAO_DISTRIB_DATA", "PREVISAO_PROD_DATA")
            VALUES (%i, %i, %i, %i, '%s', '%s', '%s');
        ''' % (row['ID_PROTOCOLO']
           , row['PREVISAO_ADMIT']
           , row['PREVISAO_DISTRIB']
           , row['PREVISAO_PROD']
           , row['PREVISAO_ADMIT_DATA']
           , row['PREVISAO_DISTRIB_DATA']
           , row['PREVISAO_PROD_DATA']
        )
    cursor.execute(qins)

# Close the cursor and connection
cursor.close()
conn.close()
