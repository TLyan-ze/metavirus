import sys
import numpy as np
import pandas as pd
from pandas import DataFrame
import openpyxl


merge1 = sys.argv[1]#1是CAT的结果信息
merge2 = sys.argv[2]#2是rdrp的结果信息
merge3 = sys.argv[3]#3是checkv的结果信息
merge4 = sys.argv[4]#4是输出文件
merge5 = sys.argv[5]#5是样本名
merge6 = sys.argv[6]#6是blastn的结果信息
merge7 = sys.argv[7]#7是sofa的深度结果信息
merge8 = sys.argv[8]#8是共线性信息
merge9 = sys.argv[9]#9是共线性覆盖度
merge10 = sys.argv[10] #10是ICTV的宿主信息
merge11 = sys.argv[11] #11是rdrp对应科信息

f1=open(merge1)
f2=open(merge2)
f3=open(merge3)

f6=open(merge6)
f7=open(merge7)

f8=open(merge8)
f9=open(merge9)

###RdRp数据框
result_dict_cat = {}
result_dict_rdrp = {}
result_dict_checkv = {}
result_dict_blastn = {}
result_dict_coverage = {}
result_dict_gong_info = {}
result_dict_gong_coverage = {}

cat_1=[]
cat_2=[]
cat_3=[]
cat_4=[]
cat_5=[]
cat_6=[]

rdrp_1=[]
rdrp_2=[]
rdrp_3=[]
rdrp_4=[]
rdrp_5=[]

checkv_1=[]
checkv_2=[]
checkv_3=[]
checkv_4=[]
checkv_5=[]#置信度


blastn_1=[]
blastn_2=[]
blastn_3=[]
blastn_4=[]
blastn_5=[]
blastn_6=[]
blastn_7=[]
blastn_8=[]
blastn_9=[]#

coverage_1=[]
coverage_2=[]
coverage_3=[]

gong_info_1=[]
gong_info_2=[]
gong_info_3=[]

gong_coverage_1=[]
gong_coverage_2=[]
##########rdrp###########################
#for line in f2.readlines():
#	line =line.strip()
        #print(line)
#	lines = line.split('\t')
	#print(lines)
#	rdrp_1.append(lines[0])
#	rdrp_2.append(lines[3])
#	rdrp_3.append(lines[4])
#	rdrp_4.append(lines[11])
#	rdrp_5.append(lines[2])
        #print(row)
#result_dict_rdrp['contig_rdrp_id'] = rdrp_1
#result_dict_rdrp['contig_identities'] = rdrp_2
#result_dict_rdrp['contig_rdrp_length'] =rdrp_3
#result_dict_rdrp['rdrp_blastp_E-value'] = rdrp_4
#result_dict_rdrp['rdrp_blastp_acc'] = rdrp_5
#print(len(contig_1))
#print(len(ident_1))
#datafa_rdrp = DataFrame(result_dict_rdrp)
#print(datafa1)
#######################################CAT
for line in f1.readlines():
	row = []  # 记录每一行
	line =line.strip()
	lines = line.split(' ')
	cat_1.append(lines[0])
	#print(lines[0])
	cat_2.append(lines[1])
	cat_3.append(lines[2])#order
	cat_4.append(lines[3])#family
	cat_5.append(lines[4])#genes
	cat_6.append(merge5)#samples
	#result_dict2[lines[0].split(' ', 1)[0]] = lines[0].split(' ', 1)[1]

result_dict_cat['contig_CAT_id'] = cat_1
result_dict_cat['contig_type'] = cat_2
result_dict_cat['contig_Order'] = cat_3
result_dict_cat['contig_Family'] = cat_4
result_dict_cat['contig_Genus'] = cat_5
result_dict_cat['sample_name']= cat_6
datafa_cat = DataFrame(result_dict_cat)

#print(datafa2)
###################rdrp
df_info= pd.DataFrame(pd.read_csv(merge2,header =None, sep= '\t', low_memory = False))
df_info2 = pd.DataFrame(pd.read_csv(merge11,sep= '\t',low_memory = False))
#print (df_info2)
#加上科的信息
datafa_samples2 = pd.merge(df_info,df_info2,left_on=2,right_on='RdRp',how='left')
#with pd.ExcelWriter('D:/work/wuhan/temp_0531.xlsx') as writer:
#    datafa_samples2.to_excel(writer, sheet_name='samples_info', index=False)
#取出比对结果为最优解
nor= datafa_samples2.shape[0]
dict1={}
dict2= {}
for i in range(0,nor):
    if datafa_samples2.loc[i,0] not in dict1.keys():
        dict1[datafa_samples2.loc[i,0]] = datafa_samples2.loc[i,12]
        dict2[datafa_samples2.loc[i,0]] = i
    elif datafa_samples2.loc[i,0] in dict1.keys():
        if dict1[datafa_samples2.loc[i,0]] >= datafa_samples2.loc[i,12]:
            print(1)
        elif dict1[datafa_samples2.loc[i,0]] < datafa_samples2.loc[i,12]:
            dict2[datafa_samples2.loc[i,0]] = i

result_dict_rdrp={}
rdrp_1=[]
rdrp_2=[]
rdrp_3=[]
rdrp_4=[]
rdrp_5=[]
rdrp_6=[]
rdrp_7=[]
for key,value in dict2.items():
    rdrp_1.append(datafa_samples2.loc[value,0])
    rdrp_2.append(datafa_samples2.loc[value,2])
    rdrp_3.append(datafa_samples2.loc[value,3])
    rdrp_4.append(datafa_samples2.loc[value,4])
    rdrp_5.append(datafa_samples2.loc[value,11])
    rdrp_6.append(datafa_samples2.loc[value,12])
    rdrp_7.append(datafa_samples2.loc[value,'family'])
result_dict_rdrp['contig_rdrp_id'] = rdrp_1
result_dict_rdrp['contig_rdrp_acc'] = rdrp_2
result_dict_rdrp['contig_identities'] = rdrp_3
result_dict_rdrp['contig_rdrp_length'] =rdrp_4
result_dict_rdrp['rdrp_blastp_E-value'] = rdrp_5
result_dict_rdrp['rdrp_blastp_sorce'] = rdrp_6
result_dict_rdrp['rdrp_blastp_family'] =rdrp_7
datafa_rdrp = DataFrame(result_dict_rdrp)
#print (datafa_rdrp)
#取出含有科
datafa_samples2['family'] = datafa_samples2['family'].replace(np.nan, 0)
csv_contig_family = datafa_samples2[~(datafa_samples2['family'].isin([0]))]
csv_contig_new = csv_contig_family.drop_duplicates(subset=[0])
datafa_rdrp_new = pd.merge(datafa_rdrp,csv_contig_new,left_on='contig_rdrp_id',right_on=0,how='left')
print (datafa_rdrp_new['contig_rdrp_id'])

########################checkv
for line in f3.readlines():
	row = []  # 记录每一行
	line =line.strip()
	lines = line.split(' ')
	checkv_1.append(lines[0])
	#print(lines[0])
	checkv_2.append(lines[1])
	checkv_3.append(lines[2])
	checkv_4.append(lines[3])
	checkv_5.append(lines[4])
        #result_dict2[lines[0].split(' ', 1)[0]] = lines[0].split(' ', 1)[1]

result_dict_checkv['contig_checkv_id'] = checkv_1
result_dict_checkv['contig_length'] = checkv_2
result_dict_checkv['checkv_quakity'] = checkv_3
result_dict_checkv['miuvig_quality'] = checkv_4
result_dict_checkv['miuvig_completeness'] = checkv_5
datafa_checkv = DataFrame(result_dict_checkv)

#######################blastn
for line in f6.readlines():
	line =line.strip()
	#print(line)
	lines = line.split('\t')
	#print(lines)
	blastn_1.append(lines[0])
	blastn_2.append(lines[1])
	blastn_3.append(lines[3])
	blastn_4.append(lines[4])
	blastn_5.append(lines[5])
	blastn_6.append(lines[9])
	blastn_7.append(lines[10])
	blastn_8.append(lines[11])
	blastn_9.append(lines[12])
#print(row)
result_dict_blastn['contig_blastn_id'] = blastn_1
result_dict_blastn['blastn_accesion'] = blastn_2
result_dict_blastn['blastn_length'] =blastn_3
result_dict_blastn['blastn_pident'] = blastn_4
result_dict_blastn['blastn_evalue'] = blastn_5
result_dict_blastn['Subject_Scientific_Name'] = blastn_6
result_dict_blastn['Subject_Common_Name'] = blastn_7
result_dict_blastn['Subject_Blast_Name'] = blastn_8
result_dict_blastn['Subject_Super_Kingdom'] = blastn_9
#print(len(contig_1))
#print(len(ident_1))
datafa_blastn= DataFrame(result_dict_blastn)
nor=datafa_blastn.shape[0]
print(nor)
datafa_blastn['blastn_filt']=''
for i in range(0, nor):
	if datafa_blastn.at[i,'Subject_Super_Kingdom'] == 'Bacteria':
		if eval(datafa_blastn.at[i,'blastn_pident']) >= 60 :
			datafa_blastn.at[i,'blastn_filt'] = 'filt_Bacteria'
		else:
			datafa_blastn.at[i,'blastn_filt'] = 'NOfilt'
	elif datafa_blastn.at[i,'Subject_Super_Kingdom'] == 'Eukaryota':
		if eval(datafa_blastn.at[i,'blastn_pident']) >= 60 :
			datafa_blastn.at[i,'blastn_filt'] = 'filt_Eukaryota'
		else:
			datafa_blastn.at[i,'blastn_filt'] = 'NOfilt'
	else:
		datafa_blastn.at[i,'blastn_filt'] = 'NOfilt'

print(datafa_blastn)
#####################################覆盖度
for line in f7.readlines():
	row = []  # 记录每一行
	line =line.strip()
	lines = line.split(' ')
	coverage_1.append(lines[0])
	#print(lines[0])
	coverage_2.append(lines[1])
	coverage_3.append(lines[2])
	#result_dict2[lines[0].split(' ', 1)[0]] = lines[0].split(' ', 1)[1]

result_dict_coverage['contig_coverage_id'] = coverage_1
result_dict_coverage['contig_Percentage'] = coverage_2
result_dict_coverage['contig_Depth'] = coverage_3
datafa_coverage = DataFrame(result_dict_coverage)


#####################################近缘物种信息
for line in f8.readlines():
	line =line.strip()
	lines = line.split('\t')
	gong_info_1.append(lines[0])
	#print(lines[0])
	gong_info_2.append(lines[1])
	gong_info_3.append(lines[2])
	#result_dict2[lines[0].split(' ', 1)[0]] = lines[0].split(' ', 1)[1]

result_dict_gong_info['cspecies_id'] = gong_info_1
result_dict_gong_info['cspecies_info'] = gong_info_2
result_dict_gong_info['cspecies_length'] = gong_info_3
datafa_gong_info = DataFrame(result_dict_gong_info)
######################################共线性coverage
for line in f9.readlines():
	line =line.strip()
	lines = line.split('\t')
	gong_coverage_1.append(lines[0])
	#print(lines[0])
	gong_coverage_2.append(lines[1])
        #result_dict2[lines[0].split(' ', 1)[0]] = lines[0].split(' ', 1)[1]

result_dict_gong_coverage['contig_species_id'] = gong_coverage_1
result_dict_gong_coverage['cspecies_coverage'] = gong_coverage_2
datafa_gong_coverage = DataFrame(result_dict_gong_coverage)

######################合并
datafa1 = pd.merge(datafa_cat,datafa_checkv,left_on='contig_CAT_id',right_on='contig_checkv_id',how='left')
datafa2 = pd.merge(datafa1,datafa_rdrp_new,left_on='contig_CAT_id',right_on='contig_rdrp_id',how='left')
datafa3 = pd.merge(datafa2,datafa_blastn,left_on='contig_CAT_id',right_on='contig_blastn_id',how='left')
datafa4 = pd.merge(datafa3,datafa_coverage,left_on='contig_CAT_id',right_on='contig_coverage_id',how='left')
datafa5 = pd.merge(datafa4,datafa_gong_info,left_on='contig_CAT_id',right_on='cspecies_id',how='left')
df_contig_info = pd.merge(datafa5,datafa_gong_coverage,left_on='contig_CAT_id',right_on='contig_species_id',how='left')


#添加自然宿主信息
#df_Host = pd.DataFrame(pd.read_excel(merge10,skiprows=0))
#df_Host = df_Host.drop_duplicates(subset=['Family'])
#datafa_Host=df_Host[['Family','Host_Source']]
#datafa7 = pd.merge(datafa6,datafa_Host,left_on='contig_Family',right_on='Family',how='left')
#datafa7.to_csv(merge4,index=False)

df_Host = pd.DataFrame(pd.read_excel(merge10,skiprows=0))
df_Host = df_Host.where(df_Host.notnull(), 'N')
df_Host1 = df_Host[~(df_Host['Genome composition'].isin(['dsDNA','ssDNA','ssDNA (+)','ssDNA (-)','ssDNA (+/-)','dsDNA-RT','N']))]

def Host(df,b):
    Order_Host = {}
    for i in df.index:
        if str(df.loc[i, b]) not in Order_Host.keys():
            Order_Host[str(df.loc[i, b])] = str(df.loc[i, 'Host Source'])
        elif str(df.loc[i, b]) in Order_Host.keys():
            list = str(df.loc[i, 'Host Source']).split(',')
            #print(list)
            for i2 in list:
                if i2 not in Order_Host[str(df.loc[i, b])]:
                    Order_Host[str(df.loc[i, b])] = Order_Host[str(df.loc[i, b])] + ',' + str(i2)
    return Order_Host
d1 = Host(df_Host1,'Order')
d2 = Host(df_Host1,'Family')
d3 = Host(df_Host1,'Genus')

nor=df_contig_info.shape[0]
for i in range(0,nor):
    if df_contig_info.loc[i,'contig_Genus'] in d3.keys():
        df_contig_info.at[i,'Host Source'] =d3[df_contig_info.loc[i,'contig_Genus']]
    elif df_contig_info.loc[i,'contig_Genus'] not in d3.keys():
        if df_contig_info.loc[i, 'contig_Family'] in d2.keys():
            df_contig_info.at[i, 'Host Source'] = d2[df_contig_info.loc[i, 'contig_Family']]
        elif df_contig_info.loc[i,'contig_Family'] not in d2.keys():
            if df_contig_info.loc[i, 'contig_Order'] in d1.keys():
                df_contig_info.at[i, 'Host Source'] = d1[df_contig_info.loc[i, 'contig_Order']]

df_contig_info.to_csv(merge4,index=False)
#标记深度异常

