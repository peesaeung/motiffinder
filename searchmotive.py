#obtain fasta input
fasta=open("/Users/peeranut/Documents/python/motif-test_sequences.fasta","r")
dna=fasta.readlines()
fasta.close()
dnatotal=len(dna)
dnaArray=[]
for i in range(1,dnatotal,2):
    dnaArray.append(dna[i])
dnaArraycount=len(dnaArray)

#put dna matrix profiling            
def Score(s):
    dnaArray2=[[],[],[],[],[]]
    for t in range(0,dnaArraycount):
        for u in range(s[t],s[t]+10):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[],[],[],[]]
    for h in range(0,4):
        for k in range(0,10):
            profile[h].append(0)
    for c in range(0,10):
        for r in range(0,dnaArraycount):
            if dnaArray2[r][c]=="A":
                profile[0][c]=profile[0][c]+1
            elif dnaArray2[r][c]=="T":
                profile[1][c]=profile[1][c]+1
            elif dnaArray2[r][c]=="C":
                profile[2][c]=profile[2][c]+1
            elif dnaArray2[r][c]=="G":
                profile[3][c]=profile[3][c]+1

    score=0
    for e in range(0,10):
        score=score+max(profile[0][e],profile[1][e],profile[2][e],profile[3][e])
    return score

#consensus
def conprint(s):
    dnaArray2=[[],[],[],[],[]]
    for t in range(0,5):
        for u in range(s[t],s[t]+10):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[],[],[],[]]
    for h in range(0,4):
        for k in range(0,10):
            profile[h].append(0)
    for c in range(0,10):
        for r in range(0,dnaArraycount):
            if dnaArray2[r][c]=="A":
                profile[0][c]=profile[0][c]+1
            elif dnaArray2[r][c]=="T":
                profile[1][c]=profile[1][c]+1
            elif dnaArray2[r][c]=="C":
                profile[2][c]=profile[2][c]+1
            elif dnaArray2[r][c]=="G":
                profile[3][c]=profile[3][c]+1  
                
    consensus=[]
    for f in range(0,10):
        if profile[0][f]>=profile[1][f]:
            if profile[0][f]>=profile[2][f]:
                if profile[0][f]>=profile[3][f]:
                    consensus.append("A")
                else:
                    consensus.append("G")
            elif profile[2][f]>=profile[3][f]:
                consensus.append("C")
            else:
                consensus.append("G")
        elif profile[1][f]>=profile[2][f]:
            if profile[1][f]>=profile[3][f]:
                consensus.append("T")
            else:
                    consensus.append("G")
        elif profile[2][f]>=profile[3][f]:
            consensus.append("C")
        else:
            consensus.append("G")
            
    constring=""
    for g in range(0,10):
        constring=constring+consensus[g]
    return constring
    
def NextVertex(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    if i<L:
        a[i+1]=1
        #a.insert(i+1,1)
        return a,i+1
    else:
        for j in range (L,1):
            if a[j]<k:
                a[j]=a[j]+1
                return a,j
    return a,0

def Bypass(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    for j in range(i,1):
        if a[j]<k:
            a[j]=a[j]+1
            return a,j
    return a,0

def BranchAndBoundMotiveSearch(t,n,l):#txn matrix size t=number of row=5 n=total nucleotide=60 ,l=motif length=10
    s=[0,0,0,0,0,0]
    bestMotif=[]
    bestScore=0
    optimisticScore=0
    i=1
    while i>0:
        if i<t:
            optimisticScore=Score(s)+(t-i)*l
            if optimisticScore<bestScore:
                (s,i)=Bypass(s,i,t,n-l)
            else:
                (s,i)=NextVertex(s,i,t,n-l)
        else:
            if Score(s)>bestScore:
                bestScore=Score(s)
                bestMotif=s
            (s,i)=NextVertex(s,i,t,n-l)
    return bestMotif

con1=BranchAndBoundMotiveSearch(5,60,10)
scoreprint=Score(con1)
con2=conprint(con1)
