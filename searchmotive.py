#obtain fasta input
fasta=open("/Users/peeranut/Documents/python/motif-test_sequences.fasta","r")
dna=fasta.readlines()
fasta.close()
dnatotal=len(dna)
dnaArray=[]
for i in range(1,dnatotal,2):
    dnaArray.append(dna[i])
dnaArraycount=len(dnaArray)

point=[1,2,3,4,5]
Score2=0
optimisticScore=0

bestMotif=[]
    
def NextVertex(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    if i<L:
        a[i+1]=1
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

def BranchAndBoundMotiveSearch(DNA,t,n,l):#txn matrix size t=number of row=5 n=total nucleotide=60 ,l=motif length=10
    s=[0,0,0,0,0]
    bestScore=0
    
    i=1
    while i>0:
        if i<t:
            optimisticScore=Score2+(t-i)*l
            if optimisticScore<bestScore:
                (s,i)=Bypass(s,i,t,n-l+1)
            else:
                (s,i)=NextVertex(s,i,t,n-l+1)
        else:
            if Score1>bestScore:
                bestScore=Score2
                bestMotif=s
            (s,i)=NextVertex(s,i,t,n-l+1)
    return bestMotif

#put dna matrix profiling            
def Score(s):
    dnaArray2=[[],[],[],[],[]]
    for t in range(0,dnaArraycount):
        for u in range(s[t],s[t]+10):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]]
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
    
#consensus
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
            
    con2=""
    for g in range(0,10):
        con2=con2+consensus[g]
    return [score,con2]

scoreprint,con1=Score(point)
    