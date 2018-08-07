import sys
import os
import time
import re
import numpy as np
import csv

n=10
copy_dir = sys.argv[1]
order = sys.argv[2]
for dir in range(1,n+1):
    if os.path.isdir(str(dir)):
        os.system('gsutil -m cp gs://'+copy_dir+str(dir)+'/*/*/*/*.log '+'gs://'+copy_dir+str(dir)+'/*/*/*/*/*.log '+ 'gs://'+copy_dir+str(dir)+'/*/*/*/*/*/*.log '+str(dir)+'/.')
    else:
        os.system('mkdir '+str(dir))
        os.system('gsutil -m cp gs://'+copy_dir+str(dir)+'/*/*/*/*.log '+'gs://'+copy_dir+str(dir)+'/*/*/*/*/*.log '+ 'gs://'+copy_dir+str(dir)+'/*/*/*/*/*/*.log '+str(dir)+'/.')

        
format = "\d+/\d+/\d+ \d+:\d+:\d+"
tr1 = "docker run"
for dir in range(1,n+1):
    dir_name = dir
    files= os.listdir(str(dir_name))
    w=open(str(dir_name)+"/"+"result.txt",'w+')
    for file in (files):
        if ("stderr" not in file) and ("result" not in file) and ("out" not in file):
             if "-" in file :
                 if "-0" in file :
                    # with open(str(dir_name)+"/"+file) as f:
                     f=open(str(dir_name)+"/"+file)
                     file = file.replace('-0','')
                 else:
                     continue 
             else :
                 f=open(str(dir_name)+"/"+file)
             #print (file)    
             f1 = f.readlines()
             l1 = f1[0]
             l2 = f1[-1]
             start=re.findall(format,l1)
             end=re.findall(format,l2)
             i=0
             while 1:
                 i=i+1
                 if tr1 in f1[i-1]:
                     call_start = re.findall(format,f1[i-1])
                     break
                 if i == len(f1):
                     break
             if i != len(f1):
                call_end = re.findall(format,f1[i])
                start1=time.strptime(start[0], "%Y/%m/%d %H:%M:%S")
                end1=time.strptime(end[0], "%Y/%m/%d %H:%M:%S")
                call_start1=time.strptime(call_start[0], "%Y/%m/%d %H:%M:%S")
                call_end1=time.strptime(call_end[0], "%Y/%m/%d %H:%M:%S")
                total = time.mktime(end1)-time.mktime(start1)
                call = time.mktime(call_end1)-time.mktime(call_start1)
                sys = total - call
                w.write(file+" total time :"+str(int(total))+"s"+'\n')
                w.write(file+" total call time :"+ str(int(call))+"s"+'\n')  
                w.write(file+" total system time :"+ str(int(sys))+"s"+'\n')
             else:
                print (str(dir))
                print (file)
             f.close()
    w.close()




#order = sys.argv[1]
#files= os.listdir(str(dir))
format = "\d+"
#x = PrettyTable()
#x.title = "Time Comparasion"
#x.field_names = ["function","time","636M 697M","1.43G 1.55G","1.77G 1.89G","237M 258","600M 663M","354M,397M","875M 778M","796M 868M","2.5G 2.65G","1.16G 1.23G"]
csv1 = open("total.csv","w+")
w1=csv.writer(csv1)
csv2 = open("call.csv","w+")
w2=csv.writer(csv2)
csv3 = open("system.csv","w+")
w3=csv.writer(csv3)
w1.writerow(["function","0.65","1.49","1.83","0.24","0.62","0.37","0.81","0.81","2.58","1.19"])
w2.writerow(["function","0.65","1.49","1.83","0.24","0.62","0.37","0.81","0.81","2.58","1.19"])
w3.writerow(["function","0.65","1.49","1.83","0.24","0.62","0.37","0.81","0.81","2.58","1.19"])
with open(order) as f:
    lines = f.readlines()
    n=len(lines)-1
    for i in range(1,n+1):
        name = lines[i-1].strip()
        #print(name)
        total_time = [name]
        call_time = [name]
        system_time = [name]
        #for file in files:
        #if "result" in file:
        for dir in range(1,11):
              with open (str(dir)+"/"+"result.txt") as f1:
                       l=f1.readlines()
                       j=-1
                       while 1:
                            j=j+1
                            if name in l[j]:
                               total = re.findall(format,l[j])
                               #print(total)
                               call = re.findall(format,l[j+1])
                               sys = re.findall(format,l[j+2])
                               total_time.append(total[0])
                               call_time.append(call[0])
                               system_time.append(sys[0])
                               break
        w1.writerow(total_time)
        w2.writerow(call_time)
        w3.writerow(system_time)
        print(name)
        #x.add_row([name,"total_time",total_time[0],total_time[1],total_time[2],total_time[3],total_time[4],total_time[5],total_time[6],total_time[7],total_time[8],total_time[9]])
        #x.add_row([" ","call_time",call_time[0],call_time[1],call_time[2],call_time[3],call_time[4],call_time[5],call_time[6],call_time[7],call_time[8],call_time[9]])
        #x.add_row([" ","system_time",system_time[0],system_time[1],system_time[2],system_time[3],system_time[4],system_time[5],system_time[6],system_time[7],system_time[8],system_time[9]])

             
csv1.close()
csv2.close()
csv3.close()
#print(x)
#data = x.get_string()

