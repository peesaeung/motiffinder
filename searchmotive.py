#obtain fasta input
fasta=open("/Users/peeranut/Documents/python/motif-test_sequences.fasta","r")
#fasta=open("/Users/peeranut/Documents/python/test.fasta","r")
dna=fasta.readlines()
fasta.close()
dnatotal=len(dna)
dnaArray=[x for x in dna if ">" not in x ]
dnaArraycount=len(dnaArray)
dnacount2=len(dnaArray[0])-1
mlength=10
#put dna matrix profiling 

#consensus
def conprint(s,k):
    dnaArray2=[]
    for i in range(0,dnaArraycount):
        dnaArray2.append([])
    for t in range(0,dnaArraycount):
        for u in range(s[t],s[t]+k):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[],[],[],[]]
    for h in range(0,4):
        for j in range(0,k):
            profile[h].append(0)
    for c in range(0,k):
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
    for f in range(0,mlength):
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
            
#    constring=""
 #   for g in range(0,6):
  #      constring=constring+consensus[g]
    print(consensus) 
    return consensus
    
def Score(s,k):#k=motif length,w=number of rows
    dnaArray2=[]
    for i in range(0,dnaArraycount):
        dnaArray2.append([])
    for t in range(0,dnaArraycount):
        for u in range(s[t],s[t]+k):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[],[],[],[]]
    for h in range(0,4):
        for j in range(0,k):
            profile[h].append(0)
    for c in range(0,k):
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
    for e in range(0,k):
        score=score+max(profile[0][e],profile[1][e],profile[2][e],profile[3][e])
    return score

def NextVertex(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    if i<L:
        a[i]=0
        return (a,i+1)
    else:
        for j in range(L,0,-1):
            if a[j-1] < k-1:
                a[j-1]=a[j-1]+1
                return (a,j)
    return (a,0)

def Bypass(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    for j in range(i,0,-1):
        if a[j-1] < k-1:
            a[j-1] = a[j-1] +1
            return (a,j)
    return (a,0)

def BranchAndBoundMotifSearch(t,n,l):#txn matrix size t=number of row=5 n=total nucleotide=60 ,l=motif length=10
    s=[0]*t
    scorenow=0
    countLoop=0
    bestMotif=[]
    bestScore=0
    optimisticScore=0
    i=1
    while i>0:
        if i<t:
            optimisticScore=Score(s,l)+(t-i)*l
            if optimisticScore<bestScore:
                s,i=Bypass(s,i,t,n-l)
                print(countLoop,optimisticScore,bestScore,(s,i),"by")
            else:
                s,i=NextVertex(s,i,t,n-l)
                print(countLoop,optimisticScore,bestScore,(s,i),"next1")
        else:
            scorenow=Score(s,l)
            if scorenow>bestScore:
                bestScore=scorenow
                bestMotif=s.copy()
            s,i=NextVertex(s,i,t,n-l)
            print(countLoop,optimisticScore,bestScore,(s,i),"next2")
        countLoop += 1
    return bestMotif,bestScore

def GreedyMotifSearch(t,n,l):
    bestMotif=[0]*t
    s=[0]*t
    bestScore=0
    countLoop=0
    for s[0] in range(0,n-l):
        for s[1] in range(0,n-l):
            scorenow=Score(s,2)
            bestScore=Score(bestMotif,2)
            countLoop+=1
            if scorenow>bestScore:
                bestScore=scorenow
                bestMotif[0]=s[0]
                bestMotif[1]=s[1]
                print(countLoop,scorenow,bestScore,s,"next1")
    s[0]=bestMotif[0]
    s[1]=bestMotif[1]
    for i in range(2,t):
        for s[i] in range(0,n-l):
            countLoop+=1
            scorenow=Score(s,l)
            bestScore=Score(bestMotif,l)
            if scorenow>bestScore:
                bestScore=scorenow
                bestMotif[i]=s[i]
                print(countLoop,scorenow,bestScore,s,"next2")
        s[i]=bestMotif[i]
    return bestMotif,bestScore

#best,scoreprint=BranchAndBoundMotifSearch(dnaArraycount,dnacount2,mlength)
#con=conprint(best)

scoreprint2=Score([23,0,42,17,19],mlength)
con2=conprint([23,0,42,17,19],mlength)

best2,scoreprint3=GreedyMotifSearch(dnaArraycount,dnacount2,mlength)
con3=conprint([0,28,8,13,25],mlength)


