import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

filename="sample_in_d_3_t_1_5_vines.txt"

disp_vine_cnt=100

dim_of_interest=[1,2]

inf_ratio=1.1

vyard={}


lines=open(filename).readlines()
max_dis=float(lines[0].strip())
# print(max_dis)

for i in range(1,len(lines)):
    line = lines[i].strip()
    if (line == ""):
        continue

    s=re.split(' ', line)
    # print(s)

    op=s[0]
    rid=int(s[1])
    b=float(s[2])
    d=float(s[3])
    dis=float(s[4])

    if dis==float('inf'):
        dis=max_dis*inf_ratio

    if op=='s':
        assert not rid in vyard
        dim=int(s[5])

        vyard[rid]={}
        vyard[rid]['dim']=dim
        vyard[rid]['rid']=rid
        vyard[rid]['start']=dis
        vyard[rid]['pts']=[(b,d,dis)]
    elif op=='c' or op=='e':
        assert rid in vyard
        vyard[rid]['pts'].append((b,d,dis))
        assert not 'len' in vyard[rid]

        if op=='e':
            vyard[rid]['len']=vyard[rid]['start']-dis
            assert vyard[rid]['len'] >= 0

vyard_filtered=[]
for v in vyard.values():
    if v['dim'] in dim_of_interest:
        vyard_filtered.append(v)

print("vine count: "+str(len(vyard)))
vyard_sorted=sorted(vyard_filtered, key=lambda x: x['len'], reverse=True)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for i in range( 0,min(len(vyard_sorted),disp_vine_cnt) ):
    print("disp "+str(vyard_sorted[i]['rid'])+" "+str(vyard_sorted[i]['dim']))

    color='C'+str(i%10)
    pts=vyard_sorted[i]['pts']

    for j in range( 0,len(pts)-1 ):
        ax.plot( [pts[j][0], pts[j+1][0]], [pts[j][1], pts[j+1][1]], color, zs=[pts[j][2], pts[j+1][2]], linewidth=2)

plt.show()
Axes3D.plot()
