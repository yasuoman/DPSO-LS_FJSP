# project : DPSO+LS_FJSP
# file   : LS.py
# author:yasuoman
# datetime:2021/5/6 17:33
# software: PyCharm

"""
description：
说明：
"""
'''
局部搜索，用于在最优解附近进行局部搜索，参考书籍：柔性作业车间调度智能算法及其应用 4.3.5领域结构研究 2移动两道工序结构
'''
import Decode as Decode
import numpy as np
import copy

#找到这个甘特图上的关键路径及其关键块
def critical_path(chr,job_op_num,p_table):
    Fit,Start_time,End_time = Decode.decode(chr,job_op_num,p_table,'decode',None)

    # 找出一条关键路径,键为机器序号，值为相应的关键块
    Key_Process = {}
    # Machine_site = []
    Total_Key_Process=[]
    Current_Machine = None
    Current_Site = None
    for i in range(len(End_time)):
        #从后向前找
        Max_time = np.max(End_time[i])
        if Max_time == Fit:
            E = list(End_time[i])
            # Machine_site.append(i)
            #最后的结束的工序的位置
            Maxtime_site = E.index(Max_time)
            #当前的site为结束工序的位置
            Current_Site = Maxtime_site
            #当前的mchine为i
            Current_Machine = i
            op= Decode.inverse_op_in_m(Current_Site, len(job_op_num), job_op_num)
            if Current_Machine in Key_Process:
                Key_Process[Current_Machine].append(op)
            else:
                Key_Process[Current_Machine]=[op]
            break
    #找到其中一个最后的结束工序的开始时间
    Start_Otime = Start_time[Current_Machine][Current_Site]
    while Start_Otime != 0:
        try:
            #去尝试在当前的这个机器上去找前一个紧邻的工序
            E_1 = list(End_time[Current_Machine])
            Current_Site = E_1.index(Start_Otime)
            Start_Otime = Start_time[Current_Machine][Current_Site]

        except:
            #最终返回的是多个关键块的一个列表，每个元素为一个列表
            #列表的值为在相应加工机器上的序号，值为每个关键块上的工序
            Total_Key_Process.append(Key_Process)
            #如果是找上一个site,那么一定是独立的关键块
            Key_Process={}
            # 这行代码会陷入死循环啊
            # if Current_Site % Max_O_len != 0:
            #从当前的site找上一个site，因为关键路径中，一定包含的是具有先后顺序的工件的工序
            Current_Site = Current_Site - 1
            for j in range(len(End_time)):
                if End_time[j][Current_Site] != 0:
                    Current_Machine = j
                    Start_Otime = Start_time[j][Current_Site]


        finally:
            op = Decode.inverse_op_in_m(Current_Site, len(job_op_num), job_op_num)
            if Current_Machine in Key_Process:
                Key_Process[Current_Machine].append(op)
            else:
                Key_Process[Current_Machine] = [op]

    #将最后的关键块装入
    Total_Key_Process.append(Key_Process)

    # print(Total_Key_Process)
    return Total_Key_Process
    # return CHS_O, Machine_site, CHS_M, Fit
    # 根据工序位置找染色体所在位置

#改变染色体中的OS
def change_OS(OS,Total_Key_Process):
    #处理第一个和最后一个关键块
    first_block = list(Total_Key_Process[0].values())[0]
    last_block = list(Total_Key_Process[-1].values())[0]
    #必须大于1个
    if len(first_block)>1:
        #直接取最后两个关键块
        block1 = first_block[-1]
        block2 = first_block[-2]
        #如果属于不同的工件才换
        if block1[0] != block2[0]:
            count=0
            for index,os in enumerate(OS):
                if os == block1[0]:
                    count+=1
                if count == block1[1]:
                    block1_index = index
                    break
            count=0
            for index, os in enumerate(OS):
                if os == block2[0]:
                    count += 1
                if count == block2[1]:
                    block2_index = index
                    break
            temp = OS[block1_index]
            OS[block1_index]=OS[block2_index]
            OS[block2_index]=temp

    if len(last_block)>1:
        #直接取前两个关键块
        block1 = last_block[0]
        block2 = last_block[1]
        #如果属于不同的工件才换
        if block1[0] != block2[0]:
            count=0
            for index,os in enumerate(OS):
                if os == block1[0]:
                    count+=1
                if count == block1[1]:
                    block1_index = index
                    break
            count=0
            for index, os in enumerate(OS):
                if os == block2[0]:
                    count += 1
                if count == block2[1]:
                    block2_index = index
                    break
            temp = OS[block1_index]
            OS[block1_index]=OS[block2_index]
            OS[block2_index]=temp
    #处理其他的关键块
    for index in range(1,len(Total_Key_Process)-1) :
        block = list(Total_Key_Process[index].values())[0]
        # 如果长度等于2，直接交换就行了
        if len(block)==2:
            block1=block[0]
            block2=block[1]
            #如果属于不同的工件才换
            if block1[0] != block2[0]:
                count=0
                for index,os in enumerate(OS):
                    if os == block1[0]:
                        count+=1
                    if count == block1[1]:
                        block1_index = index
                        break
                count=0
                for index, os in enumerate(OS):
                    if os == block2[0]:
                        count += 1
                    if count == block2[1]:
                        block2_index = index
                        break
                temp = OS[block1_index]
                OS[block1_index]=OS[block2_index]
                OS[block2_index]=temp
        #如果长度大于2，就要交换首位四个元素
        if len(block)>2:
            first_1 = block[0]
            first_2 = block[1]

            last_1 = block[-1]
            last_2 = block[-2]

            if first_1[0]!=first_2[0]:
                count = 0
                for index, os in enumerate(OS):
                    if os == first_1[0]:
                        count += 1
                    if count == first_1[1]:
                        first_1_index = index
                        break
                count = 0
                for index, os in enumerate(OS):
                    if os == first_2[0]:
                        count += 1
                    if count == first_2[1]:
                        first_2_index = index
                        break
                temp = OS[first_1_index]
                OS[first_1_index] = OS[first_2_index]
                OS[first_2_index] = temp
            if last_1[0] != last_2[0]:
                count = 0
                for index, os in enumerate(OS):
                    if os == last_1[0]:
                        count += 1
                    if count == last_1[1]:
                        last_1_index = index
                        break
                count = 0
                for index, os in enumerate(OS):
                    if os == last_2[0]:
                        count += 1
                    if count == last_2[1]:
                        last_2_index = index
                        break
                temp = OS[last_1_index]
                OS[last_1_index] = OS[last_2_index]
                OS[last_2_index] = temp

    return OS

def LS(chr,job_op_num,p_table):
    Total_Key_Process = critical_path(chr, job_op_num, p_table)
    copy_chr = copy.deepcopy(chr)
    MS = copy_chr[:p_table.shape[0]]
    before_OS = copy_chr[p_table.shape[0]:]
    after_OS = change_OS(before_OS, Total_Key_Process)
    after_chr = np.hstack((MS, after_OS))
    before_fitness= Decode.decode(chr,job_op_num,p_table,'decode',None)[0]
    after_fitness = Decode.decode(after_chr, job_op_num, p_table, 'decode', None)[0]
    # print('before:',before_fitness)
    # print('after:',after_fitness)
    if after_fitness<before_fitness:
        return after_chr
    else:
        return chr