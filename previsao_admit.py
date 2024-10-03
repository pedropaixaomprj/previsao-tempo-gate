# Autor: Pedro Paixão
# Data: 19/09/2024

##################### PREVISÃO DE TEMPO ATÉ CONCLUSÃO DE SAT ######################

# Esse código tem como função realizar previsões de tempo até a conclusão de uma SAT.
# Engloba três tipos de cálculo diferentes: um para SATs que estão na fase de admissibilidade,
# outro para SATs na fase de admissão e finalmente um último para SATs já distribuídas.
# O código será apresentado em módulos, com cada parte sendo descrita de forma geral.
# Para maiores detalhes, recomendamos a leitura do texto completo de descrição do algoritmo.

##################################################################################
############################## IMPORTAÇÃO DE PACOTES #############################
##################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
import math
import bisect
import os

filePath=os.path.expanduser('~')+'/filesUpload/'

##################################################################################
############################### LEITURA DE ARQUIVOS ##############################
##################################################################################

# Aqui vou ter que atualizar para ler direto do banco... Ainda não consegui fazer isso pelo Spoon porque não está aceitando o código em Python.

admit = pd.read_csv(filePath+'mvw_admit.csv')
distrib = pd.read_csv(filePath+'mvw_distrib.csv')
prod = pd.read_csv(filePath+'mvw_prod.csv')
tps = pd.read_csv(filePath+'funcs_atual.csv', sep =';')
feriados = pd.read_csv(filePath+'feriados_escalas_fds.csv')
afasts = pd.read_csv(filePath+'afastamentos_funcs.csv')

# Muda tipos para "datetime" (será necessário para realizar algumas transformações futuras)
afasts.loc[:, 'INICIO'] = pd.to_datetime(afasts.loc[:, 'INICIO'])
afasts.loc[:, 'FIM'] = pd.to_datetime(afasts.loc[:, 'FIM'])
prod['CRIACAO_IT'] = pd.to_datetime(prod['CRIACAO_IT'])
feriados['DATA'] = pd.to_datetime(feriados['DATA'], format='%d/%m/%Y')

##################################################################################
############################# DEFINE PARÂMETROS ##################################
##################################################################################

# Lista temas a serem incluídos
temas = sorted(pd.Series([item for sublist in prod['TEMAS'].str.split(',') for item in sublist]).unique())
temas = [tema for tema in temas if tema not in ['PSICOLOGIA', '###   SEM TEMA   ###']]

# Define percentil que será utilizado para obter estatística do tempo entre distribuição e criação da IT para cada tema
percentil_tempo = 0.5

# Define dias a considerar no filtro do tema
dias_a_considerar = 180

#########################################################################################################
####### CALCULA TEMPO DE ESCRITA E VELOCIDADE DA FILA PARA ACPs e NÃO-ACPs USANDO DIAS_PROD_IT ##########
#########################################################################################################

# Tempo de escrita se refere ao tempo entre distribuição e produção
# tempo_escrita_temas_nACP = {}
# tempo_escrita_temas_ACP = {}
# 
# for tema in temas:
#  
#    prod_tema_recente = prod.loc[(prod.TEMAS.str.contains(tema)) & (pd.to_datetime(prod.CRIACAO_IT) >= datetime.today() - timedelta(days=dias_a_considerar))]
#    n_funcionarios_tema = len(tps.loc[tps.TEMA == tema])
# 
#    if len(prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 10:  # Ao menos 5 ITs produzidas nos últimos 180 dias
#
#        # Cálculo usando apenas produções recentes do tema
#        tempo_escrita_temas_nACP[tema] = prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)
#        
#    else:
#
#        # Cálculo histórico, considerando o número de funcionários atuais no tema
#        tempo_escrita_temas_nACP[tema] = prod.loc[prod.TIPO_PRAZO != 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)
#
#    # Geralmente tem poucas ACPs, então usamos o tempo histórico de produção média
#    # Cálculo histórico, considerando o número de funcionários atuais no tema
#    prod_tema = prod.loc[prod.TEMAS.str.contains(tema)]
#
#    if len(prod_tema_recente.loc[prod_tema_recente.TIPO_PRAZO != 'PRAZO PROCESSUAL']) > 5:
#        tempo_escrita_temas_ACP[tema] = prod_tema.loc[prod_tema.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)
#
#    else:
#        tempo_escrita_temas_ACP[tema] = prod.loc[prod.TIPO_PRAZO == 'PRAZO PROCESSUAL'].DIAS_PROD_IT.quantile(percentil_tempo)

##################################################################################################
################ CRIANDO DICIONÁRIO DO Nº DE TPs ATIVOS EM CADA DIA POR TEMA #####################
##################################################################################################

N_tps_temas_datas = {}

for tema in temas:

    tps_tema = tps.loc[(tps.TEMA == tema) & (tps.STATUS == 1)]
    N_tps_tema = len(tps_tema)
    afasts_tema = afasts.loc[afasts.TEMA == tema]
    
    start_date = pd.to_datetime('2024-01-01')
    end_date = datetime.today() + timedelta(days=730)
    delta = timedelta(days=1)

    N_tps_tema_datas = []
    
    while start_date <= end_date:

        # Filter rows where start_date falls between INICIO and FIM
        N_tps_afastados_tema = len(afasts_tema.loc[(afasts_tema['INICIO'] <= start_date) & (start_date <= afasts_tema['FIM']), 'ID_USUARIO'].unique())
        N_tps_ativos_tema = N_tps_tema - N_tps_afastados_tema
        N_tps_tema_datas.append(N_tps_ativos_tema)
        start_date += delta
        
    N_tps_temas_datas[tema] = N_tps_tema_datas

# Transformando dicionário em DataFrame

start_date = pd.to_datetime('2024-01-01')
end_date = datetime.today() + timedelta(days=730)
delta = timedelta(days=1)

datas = [(start_date + timedelta(days=i)).strftime('%Y-%m-%d') for i in range((end_date - start_date).days + 1)]

N_tps_temas_datas_df = pd.DataFrame({'DATA': datas})

for tema in temas:
    N_tps_temas_datas_df[tema] = N_tps_temas_datas[tema]
    N_tps_temas_datas_df[tema] = N_tps_temas_datas_df[tema].astype('float')

N_tps_temas_datas_df.set_index('DATA')

##################################################################################################
############### AJUSTANDO A PRODUTIVIDADE DE COORDENADORES E TEMAS SECUNDÁRIOS ###################
##################################################################################################

# Retirar ou alterar essa seção caso dinâmica de coordenadores mude.

coords = ['DANIEL FONTANA OBERLING', 'ERIKA CANTANHEDE WUILLAUME', 'IZABELLA KRAICHETE LENTINO BARANDIER', 'LEONARDO DE SOUZA DA CONCEIÇÃO'] 

coords_escalas = {'DANIEL FONTANA OBERLING': ('VALORAÇÃO DANOS AMBIENTAIS', -0.85),
                  'ERIKA CANTANHEDE WUILLAUME': ('ENGENHARIA QUÍMICA', -0.85),
                  'IZABELLA KRAICHETE LENTINO BARANDIER': ('MOBILIDADE E TRANSPORTE', -0.85),
                  'LEONARDO DE SOUZA DA CONCEIÇÃO': ('ORÇAMENTO', -0.85)}

substs = ['ARMANDO NOGUEIRA DA GAMA LAMELA MARTINS', 'RODRIGO VALENTE SERRA', 'VIVIAN VICENTE DE ALMEIDA', 'ROMULO LUCIANO INOCENCIO']

substs_escalas = {'ARMANDO NOGUEIRA DA GAMA LAMELA MARTINS': ('VALORAÇÃO DANOS AMBIENTAIS', 0.15),
                 'RODRIGO VALENTE SERRA': ('VALORAÇÃO DANOS AMBIENTAIS', 0.15),
                 'VIVIAN VICENTE DE ALMEIDA': ('VALORAÇÃO DANOS AMBIENTAIS', 0.15),
                 'ROMULO LUCIANO INOCENCIO': ('ORÇAMENTO', 0.4)}

for coord in coords:
    
    tema_coord = tps.loc[tps.TP == coord, 'TEMA'].values[0]
    afasts_coord = afasts.loc[afasts.TP == coord]
    
    for data in datas:
        if len(afasts_coord.loc[(afasts_coord['INICIO'] <= pd.to_datetime(data)) & (pd.to_datetime(data) <= afasts_coord['FIM']), 'ID_USUARIO']) == 0: # Se não está afastado
            N_tps_temas_datas_df.loc[N_tps_temas_datas_df.DATA == data, tema_coord] += coords_escalas[coord][1]

for sub in substs:

    tema_original_sub = tps.loc[tps.TP == sub, 'TEMA'].values[0]
    tema_sec_sub = substs_escalas[sub][0]
    afasts_sub = afasts.loc[afasts.TP == sub]

    for data in datas:
        if len(afasts_sub.loc[(afasts_sub['INICIO'] <= pd.to_datetime(data)) & (pd.to_datetime(data) <= afasts_sub['FIM']), 'ID_USUARIO']) == 0: # Se não está afastado
            N_tps_temas_datas_df.loc[N_tps_temas_datas_df.DATA == data, tema_original_sub] -= substs_escalas[sub][1]
            N_tps_temas_datas_df.loc[N_tps_temas_datas_df.DATA == data, tema_sec_sub] += substs_escalas[sub][1]


##################################################################################################
############### AJUSTANDO A PRODUTIVIDADE EM FINAIS DE SEMANA E FERIADOS ###################
##################################################################################################

# Filtra linhas em 'feriados' nas quais 'TIPO_DIA' is 'Feriado' or 'Escala'

feriados_escalas = pd.read_csv(filePath+'feriados_escalas_fds.csv')
feriados_escalas['DATA'] = pd.to_datetime(feriados_escalas['DATA'], format = '%d/%m/%Y').dt.strftime("%Y-%m-%d")
# feriados_escalas = feriados_escalas.loc[(feriados_escalas.TIPO_DIA.isin(['Feriado', 'Escala', 'Fim de semana'])) & (feriados_escalas.DATA.isin(N_tps_temas_datas_df.DATA))]
feriados_escalas = feriados_escalas.loc[(feriados_escalas.TIPO_DIA.isin(['Feriado', 'Escala'])) & (feriados_escalas.DATA.isin(N_tps_temas_datas_df.DATA))]

for tema in temas:
    for data in feriados_escalas['DATA'].values:
            N_tps_temas_datas_df.loc[N_tps_temas_datas_df.DATA == data, tema] *= feriados_escalas.loc[feriados_escalas.DATA == data]['FATOR_PROD'].values

##################################################################################################
##################### ALTERA PLANILHA DE ADMITIDOS DA GESTÃO DO ACERVO ###########################
################################## SEPARA ACPs E NÃO-ACPs ########################################
##################################################################################################

dic_filas_com_previsao = {}

for tema in temas:
    
    # Pega apenas SEIs no tema e limpa duplicatas
    admitt = admit[admit['TEMAS'].str.contains(tema)].drop_duplicates('SEI')
    
    # Obtém posição no tema
    admitt['POS_TEMA'] = np.arange(1,len(admitt)+1)

    # Adiciona um identificador para ACP e não-ACP
    admitt['ACP'] = admitt['TIPO_PRAZO'].apply(lambda x: 'S' if x == 'PRAZO PROCESSUAL' else 'N')

    # Adiciona posição em relação as ACPs e não-ACPs
    
    pos_ACP = []
    pos_nACP = []
    i_ACP = 1
    i_nACP = 1
    i_real = 1
    
    for index, row in admitt.iterrows():
        if row['ACP'] == 'N':
            pos_nACP.append(i_nACP)
            pos_ACP.append(np.nan)
            i_nACP += 1
        else:
            pos_nACP.append(np.nan)
            pos_ACP.append(i_ACP)
            i_ACP += 1

    admitt['POS_TEMA_ACP'] = pos_ACP
    admitt['POS_TEMA_nACP'] = pos_nACP

    # Para obter a posição real, as ACPs entram na frente da fila
    # O código a seguir reorganiza a fila, colocando todas as ACPs nas posiçoes da frente
    admitt.loc[admitt['POS_TEMA_ACP'].notna(), 'POS_REAL'] = admitt['POS_TEMA_ACP']

    # Determinar o número de ACPs
    k = admitt['POS_TEMA_ACP'].notna().sum()

    # Posicionando as não-ACPs
    admitt.loc[admitt['POS_TEMA_nACP'].notna(), 'POS_REAL'] = admitt['POS_TEMA_nACP'] + k

    # Ordenando o dataset pela posição real
    admitt.sort_values(by='POS_REAL', inplace=True)

    dic_filas_com_previsao[tema] = admitt

################################################################################################################
####### CALCULA TEMPO DE ESCRITA E VELOCIDADE DA FILA PARA ACPs e NÃO-ACPs USANDO NÚMERO DE SATs PRODUZIDAS ####
################################################################################################################

# Tempo de escrita se refere ao tempo entre distribuição e produção
tempo_escrita_temas_nACP = {}
tempo_escrita_temas_ACP = {}

for tema in temas:
    
    data_minima_estimacao_prod = pd.to_datetime('2024-08-01')
    data_maxima_estimacao_prod = pd.to_datetime('2024-10-01')
    dias_estimacao_prod = (data_maxima_estimacao_prod - data_minima_estimacao_prod).days

    # Não dá pra estimar condições higiênico sanitárias, porque só tem 1 TP e quase não chega produção para esse tema
    # Para fazer a previsão, vou usar o outro tema
    if tema == 'CONDIÇÕES HIGIÊNICO-SANITÁRIAS':
        tempo_escrita_temas_nACP[tema] = 60

    else:
        
        # Previsão utilizando os últimos dias (dias_a_considerar)

        # dias_estimacao_prod = 180
        
        # prod_tema_recente = prod.loc[(prod.TEMAS.str.contains(tema)) & 
        #                             (pd.to_datetime(prod.CRIACAO_IT) >= datetime.today() - timedelta(days=dias_estimacao_prod))]

        # Olha a média de funcionários ativos no tempo
        # n_funcionarios_tema_recente = np.mean(N_tps_temas_datas_df.loc[
        #                                       pd.to_datetime(N_tps_temas_datas_df.DATA) >= datetime.today() - timedelta(days=dias_estimacao_prod),
        #                                       tema])

        # n_its_produzidas_recente = len(prod_tema_recente)

        # tempo_escrita_temas_nACP[tema] = dias_estimacao_prod * n_funcionarios_tema_recente / n_its_produzidas_recente

        # Previsão atual, usando apenas a produtividade entre 01/08 e 01/10
        
        prod_tema_recente = prod.loc[(prod.TEMAS.str.contains(tema)) & 
                                     (pd.to_datetime(prod.CRIACAO_IT) >= data_minima_estimacao_prod) &
                                     (pd.to_datetime(prod.CRIACAO_IT) <= data_maxima_estimacao_prod) ]
       
        n_funcionarios_tema_recente = np.mean(N_tps_temas_datas_df.loc[
                                              (pd.to_datetime(N_tps_temas_datas_df.DATA) >= data_minima_estimacao_prod) &
                                              (pd.to_datetime(N_tps_temas_datas_df.DATA) <= data_maxima_estimacao_prod),
                                              tema])
        
        n_its_produzidas_recente = len(prod_tema_recente)
        
        tempo_escrita_temas_nACP[tema] = dias_estimacao_prod * n_funcionarios_tema_recente / n_its_produzidas_recente
           
##################################################################################################
############################## CALCULA TEMPO DE DISTRIBUIÇÃO #####################################
##################################################################################################

# Não foi observada muita diferença de tempo de produção entre ACPs e não-ACPs.
# Por isso, nesse caso, usaremos apenas as ACPs

ITs_produzidas_temas_datas = {tema: list(N_tps_temas_datas_df[tema]/tempo_escrita_temas_nACP[tema]) for tema in temas}
# ITs_produzidas_temas_datas_ACP = {tema: list(N_tps_temas_datas_df[tema]/tempo_escrita_temas_ACP[tema]) for tema in temas}

# É um somatório acumulado. Em termos contínuos, seria o equivalente de uma integral de 0 (data atual) até o turning point.
# O resultado disso será uma sequência que dirá quantas ITs foram produzidas até o fim de cada dia.
# Por exemplo, para um tema no qual se produzem 0.6 ITs por dia, isso seria
# [0.6, 1.2, 1.8, 2.4, 3.0, 3.6, 4.2...]
ITs_produzidas_temas_datas_acc = {tema: np.cumsum(ITs_produzidas_temas_datas[tema]) for tema in temas}

# Calcula "turning point" (isto é, o dia no qual cada IT terminou de ser produzida).
# Para o exemplo anterior, seriam os dias 2, 4, 5, 7 e assim por diante

turning_point_temas_datas = {}

for tema in temas:

    turning_point_tema_datas = []
    
    for i in range (1, len(ITs_produzidas_temas_datas_acc[tema])+1):
        index = bisect.bisect_right(ITs_produzidas_temas_datas_acc[tema], i)
        turning_point_tema_datas.append(index) # Talvez aqui colocar 1 dia a mais de garantia?
        
        # primeira posição para a qual o elemento de ITs_produzidas_temas_datas_acc[tema] é maior ou igual a i
        
    turning_point_temas_datas[tema] = turning_point_tema_datas

# Cria dicionário contendo dataframes por tema, com a previsão

##################################################################################################
########################      ESTIMA CHEGADA DE ACP's QUE FURAM A FILA ###########################
##################################################################################################

# Essa seção tem como objetivo obter estatísticas sobre a chegada de ACPs nos diferentes temas ao longo do tempo.
# Em particular, queremos saber qual é a porcentagem de ACPs dentre o total de procedimentos. Esse parâmetro
# será usado como um "multiplicador" ao fim do processo para adicionar ao tempo previsto. Por exemplo, se temos
# 10% de ACPs em um certo tema, a previsão final de tempo para SATs desse tema será acrescida de 10%.
# Note que alguns temas historicamente não tem ACPs. Isso pode ser ocasionado pela baixa quantidade de SATs no tema.
# Por isso, definimos uma taxa mínima de 5%.

# Por enquanto, não vou usar fator de correcao
# threshold_acp_min = 0
# threshold_acp_max = 1

#  # Calcula proporção de ACPs dentre o total no tema
#  dic_prop_ACPs_tema = {tema: min(threshold_acp_max,
#                             max(threshold_acp_min,
#                             len(prod.loc[(prod.TIPO_PRAZO == 'PRAZO PROCESSUAL') & (prod.TEMAS.str.contains(tema)), :]) / len(prod.loc[prod.TEMAS.str.contains(tema), :])))
#     for tema in temas
# }

# Calcula taxa de diária chegada de SATs admitidas no tema, utilizando os últimos 365 dias

# dias_a_considerar_taxa = 365
# data_min = datetime.today() - timedelta(days=dias_a_considerar_taxa)
# dic_taxa_chegada_ACPs_temas = {}

# for tema in temas:
    
#    n_admit_recente = len(admit.loc[(pd.to_datetime(admit.ENTRADA) >= data_min)  & (prod.TEMAS.str.contains(tema))])
#    n_distrib_recente = len(distrib.loc[(pd.to_datetime(distrib.ENTRADA) >= data_min) & (prod.TEMAS.str.contains(tema))])
#    n_prod_recente = len(prod.loc[(pd.to_datetime(prod.ENTRADA) >= data_min) & (prod.TEMAS.str.contains(tema))])
#    taxa_chegada_tema = (n_admit_recente + n_distrib_recente + n_prod_recente)/dias_a_considerar_taxa

#    # Dicionário com a taxa de chegadas de ACPs por tema
#    dic_taxa_chegada_ACPs_temas[tema] = dic_prop_ACPs_tema[tema] * taxa_chegada_tema

# Se quiser ver o DataFrame resultante
# taxa_chegada_ACPs_temas = pd.DataFrame(list(dic_taxa_chegada_ACPs_temas.items()), columns=['Tema', 'Taxa diária de chegada de ACPs'])

# Define thresholds para fator de correção

# threshold_fator_correcao_min = 0
# threshold_fator_correcao_max = 1

# Define fator de correção
# fator_correcao_temas = {tema:
#                        min(threshold_fator_correcao_max,
#                        max(threshold_fator_correcao_min,
#                        dic_taxa_chegada_ACPs_temas[tema] * tempo_escrita_temas_nACP[tema] / np.mean(N_tps_temas_datas[tema])))
#                        for tema in temas}


##################################################################################################
##################################     PREVISÃO FINAL ############################################
##################################################################################################

for tema in temas:

    # Note que a aqui dic_filas_previsao já precisa estar ordenada por "POS_REAL"
    dic_filas_com_previsao[tema]['PREVISAO_FILA'] = turning_point_temas_datas[tema][0:len(dic_filas_com_previsao[tema])]
    
    # Adicionando fator de entrada de ACPs (se quiser, retirar essa linha em específico depois)
    # dic_filas_com_previsao[tema]['PREVISAO_FILA'] = np.ceil(
    #                                                 (1 + fator_correcao_temas[tema]) *
    #                                                 dic_filas_com_previsao[tema]['PREVISAO_FILA'])

    dic_filas_com_previsao[tema]['PREVISAO_PROD'] = (tempo_escrita_temas_nACP[tema] +
                                                     dic_filas_com_previsao[tema]['PREVISAO_FILA']).astype(int)

previsao_final = pd.concat(dic_filas_com_previsao.values(), ignore_index=True)

############### Cria intervalo de segurança usando 20% a mais e a menos de produtividade
# Intervalo de tamanho no mínimo 10, mínimo de tamanho 0, máximo de tamanho 10

previsao_final['PREVISAO_FILA_MIN'] = np.maximum(0, np.minimum(previsao_final['PREVISAO_FILA'] - 5,
                                                                np.floor(0.8 * previsao_final['PREVISAO_FILA'])).astype(int))
 
previsao_final['PREVISAO_FILA_MAX'] = np.maximum(10, np.maximum(previsao_final['PREVISAO_FILA'] + 5,
                                                                 np.ceil(1.2 * previsao_final['PREVISAO_FILA'])).astype(int))

############### Cria intervalo de confiança usando 20% a mais ou a menos de produtividade

# previsao_final['PREVISAO_PROD_MIN'] = (tempo_escrita_temas_nACP[tema] + previsao_final['PREVISAO_FILA_MIN']).astype(int)

# previsao_final['PREVISAO_PROD_MAX'] = (tempo_escrita_temas_nACP[tema] + previsao_final['PREVISAO_FILA_MAX']).astype(int)


#################### CORREÇÃO FUTURA A SER ADICIONADA ########################

    # Se tem TP livre, a previsão do dia inicial tem que ser 0. Com isso, temos que fazer PREVISAO_FILA = PREVISAO_FILA - PREVISAO_FILA[0]

##################################################################################################
############# AJUSTE PARA PREVISÕES MULTITEMA - USA O PIOR DOS CASOS #############################
##################################################################################################

# For the duplicated rows, find the max of 'PREVISAO_FILA' and 'PREVISAO_PROD' for each group
max_vals = previsao_final.groupby(['SEI', 'ENTRADA', 'TEMAS'], as_index=False)[['PREVISAO_FILA', 'PREVISAO_PROD']].max()

# Merge the max values back into the original DataFrame
previsao_final[['PREVISAO_FILA', 'PREVISAO_PROD']] = previsao_final[['SEI', 'ENTRADA', 'TEMAS']].merge(max_vals, on=['SEI', 'ENTRADA', 'TEMAS'], how='left')[['PREVISAO_FILA', 'PREVISAO_PROD']]


##################################################################################################
#################################### ACRESCENTA DATAS ############################################
##################################################################################################

previsao_final['PREVISAO_DISTRIB_DATA'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_FILA'], unit='D')
previsao_final['PREVISAO_DISTRIB_DATA_MIN'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_FILA_MIN'], unit='D')
previsao_final['PREVISAO_DISTRIB_DATA_MAX'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_FILA_MAX'], unit='D')
previsao_final['PREVISAO_PROD_DATA'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_PROD'], unit='D')
# previsao_final['PREVISAO_PROD_DATA_MIN'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_PROD_MIN'], unit='D')
# previsao_final['PREVISAO_PROD_DATA_MAX'] = pd.to_datetime(datetime.today()) + pd.to_timedelta(previsao_final['PREVISAO_PROD_MAX'], unit='D')

# Função para transformar finais de semana em segundas

def proximo_dia_util(date):
    if pd.isna(date):  # Se é nan
        return date
    if date.weekday() >= 5:  # Se é sábado (5) ou domingo (6)
        return date + pd.offsets.BDay(1)
    return date  # Se já é um dia útil

##################################################################################################
################################ CONCATENA DE NOVO E EXPORTA #####################################
##################################################################################################

# Aplica a função e transforma em string
previsao_final['PREVISAO_DISTRIB_DATA'] = previsao_final['PREVISAO_DISTRIB_DATA'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d')
previsao_final['PREVISAO_DISTRIB_DATA_MIN'] = previsao_final['PREVISAO_DISTRIB_DATA_MIN'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d')
previsao_final['PREVISAO_DISTRIB_DATA_MAX'] = previsao_final['PREVISAO_DISTRIB_DATA_MAX'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d')
previsao_final['PREVISAO_PROD_DATA'] = previsao_final['PREVISAO_PROD_DATA'].apply(proximo_dia_util).dt.strftime('%Y-%m-%d')

previsao_final.to_csv(filePath+'previsao_fila_admitidos.csv', index=False)
