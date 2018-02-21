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
#def Score(s, DNA, k):
  #  """ 
   #     compute the consensus SCORE of a given k-mer 
    #    alignment given offsets into each DNA string.
     #       s = list of starting indices, 1-based, 0 means ignore
      ##     k = Target Motif length
    #"""
 #   score = 0
  #  for i in range(k):
         #loop over string positions
    #    cnt = dict(zip("ATCG",(0,0,0,0)))
     #   for j, sval in enumerate(s):
            # loop over DNA strands
       #     base = DNA[j][sval+i] 
        #    cnt[base] += 1
        #score = score+max(cnt.values())
    #return score


#consensus
def conprint(s):
    dnaArray2=[[],[],[],[]]
    for t in range(0,4):
        for u in range(s[t],s[t]+mlength):
            dnaArray2[t].append(dnaArray[t][u])

    profile=[[],[],[],[]]
    for h in range(0,4):
        for k in range(0,mlength):
            profile[h].append(0)
    for c in range(0,6):
        for r in range(0,4):
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
    
def Score2(s,k):#k=motif length,w=number of rows
    dnaArray2=[]
    for i in range(0,dnaArraycount):
        dnaArray2.append([])
    for t in range(0,dnaArraycount):
        for u in range(s[t],s[t]+k):
            dnaArray2[t].append(dnaArray[t][u])

    profile2=[[],[],[],[]]
    for h in range(0,4):
        for j in range(0,k):
            profile2[h].append(0)
    for c in range(0,k):
        for r in range(0,dnaArraycount):
            if dnaArray2[r][c]=="A":
                profile2[0][c]=profile2[0][c]+1
            elif dnaArray2[r][c]=="T":
                profile2[1][c]=profile2[1][c]+1
            elif dnaArray2[r][c]=="C":
                profile2[2][c]=profile2[2][c]+1
            elif dnaArray2[r][c]=="G":
                profile2[3][c]=profile2[3][c]+1
                
    score=0
    for e in range(0,k):
        score=score+max(profile2[0][e],profile2[1][e],profile2[2][e],profile2[3][e])
    return score

def NextVertex(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    a=[0]+a
    if i<L:
        a[i]+1
        a.pop(0)
        return (a,i+1)
    else:
        for j in range(L,0,-1):
            if a[j] < k:
                a[j]=a[j]+1
                a.pop(0)
                return (a,j)
    a.pop(0)
    return (a,0)

def Bypass(a,i,L,k): #a=array position,i=position,L=interesting line=10,k=alphabet=4
    a=[0]+a
    for j in range(i,0,-1):
        if a[j] < k:
            a[j] = a[j] +1
            a.pop(0)
            return (a,j)
    a.pop(0)
    return (a,0)

def BranchAndBoundMotifSearch2(t,n,l):#txn matrix size t=number of row=5 n=total nucleotide=60 ,l=motif length=10
    s=[]
    for i in range(0,dnaArraycount):
        s.append(0)
    countLoop=0
    bestMotif=[]
    bestScore=0
    optimisticScore=0
    i=1
    while i>0:
        if i<t:
            optimisticScore=Score2(s,l)+(t-i)*l
            if optimisticScore<bestScore:
                s,i=Bypass(s,i,t,n-l)
                print(countLoop,optimisticScore,bestScore,(s,i),"by")
            else:
                s,i=NextVertex(s,i,t,n-l)
                print(countLoop,optimisticScore,bestScore,(s,i),"next1")
        else:
            if Score2(s,l)>bestScore:
                bestScore=Score2(s,l)
                bestMotif=s
            s,i=NextVertex(s,i,t,n-l)
            print(countLoop,optimisticScore,bestScore,(s,i),"next2")
        countLoop += 1
    return bestMotif,bestScore

def GreedyMotifSearch(t,n,l):
    bestMotif=[]
    for i in range(0,dnaArraycount):
        bestMotif.append(0)
    s=[]
    for i in range(0,dnaArraycount):
        s.append(0)
    for s[0] in range(0,n-l):
        for s[1] in range(0,n-l):
            if Score2(s,l)>Score2(bestMotif,l):
                bestMotif[0]=s[0]
                bestMotif[1]=s[1]
    s[0]=bestMotif[0]
    s[1]=bestMotif[1]
    for i in range(2,t):
        for s[i] in range(0,n-l):
            if Score2(s,l)>Score2(bestMotif,l):
                bestMotif[i]=s[i]
        s[i]=bestMotif[i]
    return bestMotif

#scoreprint2=Score([12,13,14,15],dnaArray,6) 
#for i in range(0,4):
    #print(dnaArray2[i])
best,scoreprint=BranchAndBoundMotifSearch2(dnaArraycount,dnacount2,mlength)
con=conprint(best)

scoreprint2=Score2([24,1,43,18,20],mlength)
con2=conprint([24,1,43,18,20])

best2=GreedyMotifSearch(dnaArraycount,dnacount2,mlength)
scoreprint3=Score2(best2,mlength)
con3=conprint([14,1,0,0,0])

#con3=ContainedConsensusMotifSearch(dnaArray,6)
#con2=conprint([12,13,14,15])
