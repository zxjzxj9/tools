#! /usr/bin/python

from lxml import etree

# This scrpit helps to extract the dos from vasprun.xml
# Produce file:
# tdos.1 total dos for spin up
# tdos.2 total dos for spin down
# pdos should assign a range of atoms, like 1-10, 2-20
# Or 1,2,3,4,5
# pdos.1 partial dos spin up, three line s,p,d
# pdos.1 partial dos spin down, three line s,p,d



if __name__ == "__main__":
    xmldata = etree.parse("vasprun.xml")
    # Get fermi level
    tmp = xmldata.xpath("//dos/i[@name='efermi']/text()")[0]
    efermi = float(tmp)
    print("Obtain the Fermi energy, %12.6f eV" %efermi)
    print("DOS data will be aligned w.r.t to Fermi energy")

    # Get the data of Total dos, spin 1
    print("Start dumping the total dos of spin up...")
    fout = open("tdos.1","w")
    tmp = xmldata.xpath("//dos/total/array/set/set[@comment='spin 1']/r/text()")
    fout.write("%12s    %12s\n" %('Band Energy','DOS'))
    for r in tmp:
        rec = r.split()
        fout.write("%12.6f    %12.6f\n" %(float(rec[0])-efermi, float(rec[1])))
    fout.close()
    print("Dumping tdos.1 completed.")

    print("Start dumping the total dos of spin down...")
    fout = open("tdos.2","w")
    tmp = xmldata.xpath("//dos/total/array/set/set[@comment='spin 2']/r/text()")
    fout.write("%12s    %12s\n" %('Band Energy','DOS'))
    for r in tmp:
        rec = r.split()
        fout.write("%12.6f    %12.6f\n" %(float(rec[0])-efermi, float(rec[1])))
    fout.close()
    print("Dumping tdos.2 completed.")

    print("Start dumping the partial dos")
    print("Please input the atom number: e.g. 1-20 or 1,2,3,4,5,7-15")
    tmp = raw_input()
    d1 = tmp.split(",")
    atomlist = []
    for rec in d1:
        if "-" in rec:
            start = int(rec.split('-')[0])
            end = int(rec.split('-')[1])
            atomlist.extend(range(start,end+1))
        else:
            atomlist.append(int(rec))
    atomlist = set(atomlist)
    print "You have chosen atoms ", atomlist
    print("Start suming and dumping the partial dos of spin up")
    dos_s1 = []
    dos_s2 = []
    dos_p1 = []
    dos_p2 = []
    dos_d1 = []
    dos_d2 = []
    elist = []
    for atom in atomlist:
        tmp = xmldata.xpath("//partial/array/set/set[@comment='ion %d']/set[@comment='spin 1']/r/text()" %atom)
        #print("//partial/array/set/set[@comment='ion %d']/set[comment='spin 1']/r/text()" %atom)
        #print tmp
        cnt = 0
        if not elist:
            for r in tmp:
                rec = r.split()
                elist.append(float(rec[0]))
                dos_s1.append(float(rec[1]))
                dos_p1.append(float(rec[2]))
                dos_d1.append(float(rec[3]))
        else:
            dos_s1[cnt]+=float(rec[1])
            dos_p1[cnt]+=float(rec[2])
            dos_d1[cnt]+=float(rec[3])
            cnt += 1
    fout = open("pdos.1","w")
    fout.write("%12s    %12s    %12s    %12s\n"  %('Band Energy', 'S', 'P', 'D'))
    for _1,_2,_3,_4 in zip(elist,dos_s1,dos_p1,dos_d1):
        fout.write("%12.6f    %12.6f    %12.6f    %12.6f\n" %(_1-efermi,_2,_3,_4))
    fout.close()
    print("Dumping pdos.1 completed.")

    elist = []
    print("Start suming and dumping the partial dos of spin down")
    for atom in atomlist:
        tmp = xmldata.xpath("//partial/array/set/set[@comment='ion %d']/set[@comment='spin 2']/r/text()" %atom)
        #print("//partial/array/set/set[@comment='ion %d']/set[comment='spin 1']/r/text()" %atom)
        #print tmp
        cnt = 0
        if not elist:
            for r in tmp:
                rec = r.split()
                elist.append(float(rec[0]))
                dos_s2.append(float(rec[1]))
                dos_p2.append(float(rec[2]))
                dos_d2.append(float(rec[3]))
        else:
            dos_s2[cnt]+=float(rec[1])
            dos_p2[cnt]+=float(rec[2])
            dos_d2[cnt]+=float(rec[3])
            cnt += 1
    fout = open("pdos.2","w")
    fout.write("%12s    %12s    %12s    %12s\n"  %('Band Energy', 'S', 'P', 'D'))
    for _1,_2,_3,_4 in zip(elist,dos_s2,dos_p2,dos_d2):
        fout.write("%12.6f    %12.6f    %12.6f    %12.6f\n" %(_1-efermi,_2,_3,_4))
    fout.close()
    print("Dumping pdos.2 completed.")
    #xmldata.xpath("//dos/")

