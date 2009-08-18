======================== README.TXT =================================
ARTICLE INFORMATION
Journal:    Physical Review Letters
All Authors:    Kallin,Ann B. /Gonz\'alez,Iv\'an /Hastings,Matthew B.
                /Melko,Roger G. /
Title:  Valence bond and von Neumann entanglement entropy in 
        Heisenberg ladders


DEPOSIT INFORMATION
Description:
README.TXT:  This text file contains information about the
included EPAPS material.

        mfig1.eps:  This figure shows the difference (Delta S) between 
        the value of the von Neumann entropy which has been converged 
        (for large m=m_max), and the value obtained by keeping various 
        m<m_max DMRG states per block.  The value of m_max is 800 for 
        N=4, and 1200 for N=5, where N is the number of legs.  Each data 
        point is taken from an independent run of the code with fixed m, 
        which is done to converge each value of S using a large number 
        of DMRG sweeps.  In this figure, the entropy values are taken at 
        x=2N^2 (as in the paper), with ladders of length L=100 for N=4, 
        and L=4N=20 for N=5.
        
        
        mfig2.eps: This figure illustrates the convergence of the entropy 
        at x=2N^2 on a N=7 ladder with L=4N sites per leg.  It shows the 
        difference (Delta S) between the converged entropy value (obtained 
        by keeping 1800 DMRG states per block) and the value obtained with
        m states, as well as the eigenvalue truncation error.  Unlike 
        mfig1.eps, this data was obtained in one DMRG run (as is typical 
        practice for obtaining a DMRG data point), since performing multiple 
        runs to converge each m as in mfig1.eps is unpractical for such large 
        systems.  
        INSET: The DMRG run, showing the von Neumann entropy for all values 
        of x over several sweeps of the system.

Total No. of Files:     3
File Names:     README.TXT, mfig1.eps, mfig2.eps
File Types:     .TXT, .EPS
Contact Information:    
        Roger Melko
        DEPARTMENT OF PHYSICS & ASTRONOMY
        UNIVERSITY OF WATERLOO
        WATERLOO, ONTARIO, CANADA, N2L 3G1
        phone: (519) 888-4567 x38406 
        fax: (519) 746-8115
===================================================================== 
