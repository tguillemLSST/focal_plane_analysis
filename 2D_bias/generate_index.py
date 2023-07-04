
rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']

with open('index.html', 'w') as f:
    f.write("<html>\n")
    f.write("<ul>\n")
    for i in range(len(rafts)):
        for j in range(len(ccds)):
            for l in range(len(amps)) :
                f.write("<IMG align=center width=350 height=300 SRC=./mean_column_"+rafts[i]+"_"+ccds[j]+"_"+amps[l]+".png>\n")
    f.write("</ul>")
    f.write("</html>")
