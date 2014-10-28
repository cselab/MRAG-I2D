from numpy import *
import matplotlib.pyplot as plt
import sys
import subprocess


def run_reference(commonarg):
    arg = commonarg + " -fmm velocity-wim -fmm-theta 0.0 -restart 0"
    print "running reference as follows: "+arg
    output=""
    try:
        output = subprocess.check_output(arg.split())
    except subprocess.CalledProcessError, e:
        output = e.output
#    print output

def run_test(commonarg, theta, sse=True):
    arg = commonarg + " -fmm-theta "+str(theta)+ " -restart 1"
    if(sse):
        arg += " -fmm velocity -core-fmm sse"
    else:
        arg += " -fmm velocity-wim"

    print "running test as follows: "+arg
    output=""
    try:
        output = subprocess.check_output(arg.split())
    except subprocess.CalledProcessError, e:
        output = e.output
#   print output

    outlines = output.split('\n')

    L1=zeros(3)
    L2=zeros(3)
    LI=zeros(3)

    for line in outlines:
        if(line.find("L1")==0):
            L1[0]=float(line.split()[1])
            L1[1]=float(line.split()[2])
            L1[2]=float(line.split()[3])
        if(line.find("L2")==0):
            L2[0]=float(line.split()[1])
            L2[1]=float(line.split()[2])
            L2[2]=float(line.split()[3])
        if(line.find("L1")==0):
            LI[0]=float(line.split()[1])
            LI[1]=float(line.split()[2])
            LI[2]=float(line.split()[3])

    return L1,L2,LI

def compile_new(order, clean=True):
    if(clean):
        arg = "make -j 4 -C ../makefiles/ clean"
        subprocess.call(arg.split())
    arg = "make -j 4 -C ../makefiles/ if2d-tests fmm-order="+str(order)
    subprocess.call(arg.split())

def printToFile(thetas, errors,filename):
    assert(len(thetas)==len(errors))
    f = open(filename,'w')
    for i,theta in enumerate(thetas):
        f.write("%.2f %.7e %.7e %.7e\n" % (theta,errors[i,0],errors[i,1],errors[i,2]))
    f.close()

def plotResults(thetas,L1,L2,LI,SSE):
    fig=plt.figure()
    N=3
    M=1
    L1p = L1.copy()
    L2p = L2.copy()
    LIp = LI.copy()
    for i in range(0,len(L1)):
        for j in range(0,3):
            L1p[i,j] = (-30 if L1[i,j]<=0 else log10(L1[i,j]))
            L2p[i,j] = (-30 if L2[i,j]<=0 else log10(L2[i,j]))
            LIp[i,j] = (-30 if LI[i,j]<=0 else log10(LI[i,j]))

    ax = fig.add_subplot(N,M,1)
    ax.plot(thetas,L1p[:,0],'-o',lw=1.5)
    if(not SSE): ax.plot(thetas,L1p[:,1],'--k')
    ax.set_ylim(-16,2)
    ax.set_ylabel('L1 error')

    ax = fig.add_subplot(N,M,2)
    ax.plot(thetas,L2p[:,0],'-o',lw=1.5)
    if(not SSE): ax.plot(thetas,L2p[:,1],'--k')
    ax.set_ylim(-16,2)
    ax.set_ylabel('L2 error')

    ax = fig.add_subplot(N,M,3)
    ax.plot(thetas,LIp[:,0],'-o',lw=1.5)
    if(not SSE): ax.plot(thetas,LIp[:,1],'--k')
    ax.set_ylim(-16,2)
    ax.set_xlabel('theta')
    ax.set_ylabel('LI error')

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

### SOME PARAMETERS FOR THIS CHECK ###
commonarg = "../makefiles/if2d-tests -study multipole -bpd 4 -lmax 4 -uniform 1 -jump 2 -nthreads 4"
theta_list = [0., 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

if(len(sys.argv[1:])!=2):
    print 'USAGE: python '+sys.argv[0]+' < order > < withSSE >'
    sys.exit(-1)

order = int(sys.argv[1])
SSE = str2bool(sys.argv[2])
print 'RUNNING ORDER '+str(order)+' and '+("with " if SSE else "without ")+"SSE"

Nthetas=len(theta_list)
L1 = zeros((Nthetas,3))
L2 = zeros((Nthetas,3))
LI = zeros((Nthetas,3))

# COMPILE AND RUN REFERENCE
compile_new(order)
run_reference(commonarg)

# RUN FOR DIFFERENT THETAS
for i,theta in enumerate(theta_list):
    L1_out,L2_out,LI_out = run_test(commonarg,theta,SSE)
    L1[i,:] = L1_out
    L2[i,:] = L2_out
    LI[i,:] = LI_out

# STORE RESULTS IN FILE
print 'PRINTING TO FILE'
filename="errors_L1_order"+str(order)+("_SSE" if SSE else "")
printToFile(theta_list,L1,filename)
filename="errors_L2_order"+str(order)+("_SSE" if SSE else "")
printToFile(theta_list,L2,filename)
filename="errors_LI_order"+str(order)+("_SSE" if SSE else "")
printToFile(theta_list,LI,filename)

# PLOT RESULTS
print 'PLOTTING'
plotResults(theta_list,L1,L2,LI,SSE)
plt.savefig("errors_order"+str(order)+("_SSE" if SSE else "")+".pdf")
#plt.show()



