function [S,Q] = genlouvainrand(B,limit)

%GENLOUVAIN  Louvain-like community detection, specified quality function.
%   Version 1.1, January 27, 2012.
%
%   [S,Q] = GENLOUVAINRAND(B) with B a matrix implements a Louvain-like greedy
%   community detection method using the modularity/quality matrix B which
%   encodes the quality function Q obtained by summing over all elements
%   B(i,j) such that nodes i and j are placed in the same community.
%   Following Blondel et al. 2008, the algorithm proceeds in two phases
%   repeated iteratively: quality is optimized by moving one node at a
%   time until no such moves improve quality; the communities found to that
%   point are then aggregated to build a new network where each node
%   represents a community.  The output vector S encodes the obtained
%   community assignments, with S(i) identifying the community to which
%   node i has been assigned.  The output Q gives the quality of the
%   resulting partition of the network.
%
%   [S,Q] = GENLOUVAINRAND(B,limit) with B a function handle such that B(i) returns
%   the ith column of the modularity/quality matrix uses this function handle until
%   the number of groups is less than the limit specified (default 10000)
%   and then builds the B matrix corresponding to the new aggregated network 
%   in subsequent passes. Used to reduce the memory footprint for large
%   networks
%
%   Example
%         k = full(sum(A));
%         twom = sum(k); 
%         B = @(v) A(:,v) - k'*k(v)/twom;
%         [S,Q] = genlouvain(B); 
%         Q = Q/twom;
%     finds community assignments for the undirected network encoded by the
%     symmetric adjacency matrix A.  For small networks, one may obtain
%     reasonably efficient results even more simply by handling the full
%     modularity/quality matrix
%         B = A - k'*k/twom;
%     instead of the function handle.  Intended use also includes the
%     "multislice" network quality function of Mucha et al. 2010, where B
%     is built in the multislice wrapper codes listed below.
%
%   See also
%       non-randperm version:       GENLOUVAIN
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%       Multislice wrappers:        MULTIORD, MULTIORDF, MULTICAT, MULTICATF
%
%   Notes:
%     The matrix B is assumed to be both symmetric and square.  This
%     assumption is not checked here.
%
%     This version can return different results from run to run because it
%     considers nodes in pseudorandom (randperm) order.  Because of the
%     potentially large number of nearly-optimal partitions (Good et al.
%     2010), one is encouraged to investigate results of repeated
%     applications of this code.
%
%     This algorithm is only "Louvain-like" in the sense that the two
%     phases are used iteratively in the same manner as in the Louvain
%     algorithm (Blondel et al. 2008).  Because it operates on a general
%     quality/modularity matrix B, it does not include any analytical
%     formulas for quickly identifying the change in modularity from a
%     proposed move nor any improved efficiency obtained by their use.
%
%     This code has a known potential problem where accumulated subtraction
%     error might lead to an infinite loop with each pass oscillating
%     between two or more partitions yet incorrectly identifying increases
%     in quality.  We believe this problem has been ameliorated by the
%     1*eps terms in lines 168 and 238.  If you encounter this problem
%     in your use, try different multipliers and please notify Peter Mucha
%     (<a href="mailto:mucha@unc.edu">mucha@unc.edu</a>).
%
%     The output Q provides the sum over the appropriate elements of B
%     without any rescaling.  As such, the example above rescales the
%     Q output by genlouvain by 2m = sum(k) so that Q <= 1.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References:
%     Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, and
%     Etienne Lefebvre, "Fast unfolding of communities in large networks,"
%     Journal of Statistical Mechanics: Theory and Experiment, P10008
%     (2008).
%
%     Fortunato, Santo, "Community detection in graphs," Physics Reports
%     486, 75-174 (2010).
%
%     Good, Benjamin H., Yves-Alexandre de Montjoye, and Aaron Clauset,
%     "Performance of modularity maximization in practical contexts,"
%     Physical Review E 81, 046106 (2010).
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     A special thank you to Stephen Reid, whose greedy.m code was the
%     original version that has over time developed into the present code.
%
%   Citation: If you use this code, please cite as
%       Inderjit S. Jutla, Peter J. Mucha and Lucas Jeub, "A generalized Louvain method
%       for community detection implemented in MATLAB,"
%       http://netwiki.amath.unc.edu/GenLouvain (2011).
% 
% 
 

%set default for maximum size of modularity matrix
if nargin <2
    limit = 10000;
end


%initialise variables
if isa(B,'function_handle')
    n=length(B(1));
    S=(1:n)';    
    M=B;
    dtot=0;
  
else
    n = length(B);
    S = (1:n)';
    M=B;
end


%Run using function handle, if provided
while (isa(M,'function_handle'))    %loop around Blondel et al's "scond phase" with function handle
      
    y = unique(S); %unique also puts elements in ascending order
    Sb=S;  
    yb = [];
   
    clocktime=clock;
    disp(['Merging ',num2str(length(y)),' communities  ',num2str(clocktime(4:6))]);
    
    dstep=1;
    
    while (~isequal(yb,y))&&(dstep/dtot>2*eps) %This is the loop around Blondel et al's "first phase"
           
        
        
        %        Q = 0;
%        %improves performance considerably if one doesn't compute modularity
%        %for the first pass (for display purposes only)
%         P = sparse(y,1:length(y),1); %Modularity Calculation
%         for i = 1:length(M(1))
%             Q = Q + (P*M(i))'*P(:,i);
%         end
%         disp([num2str(length(unique(y))),' ',num2str(Q)])
        t=length(unique(y));
        disp(num2str(length(unique(y))))
        yb = y;
        dstep=0;
     
        
        for i = randperm(length(M(1)))    %loop over nodes in pseudorandom order     

            Mi = M(i);
          
            u = unique([y(i);y(Mi>0)]);
            
          
            dH=modchange_y(Mi,y,u);  %change in modularities calculated using y
            %implemented in C for speed (same code in matlab is slow)
           
            
            yi=find(u==y(i));
            dH(yi) = dH(yi) - Mi(i);
            
            [~, k] = max(dH);
            
            %%only move to different group if it is more optimized than
            %%staying in same group (up to error with double precision)
            if(dH(k)>(dH(yi)))
                dtot=dtot+dH(k)-dH(yi);
                dstep=dstep+dH(k)-dH(yi);
                y(i) = u(k); 
            end
            
        end
        
        disp(['change: ',num2str(dstep),' total: ',num2str(dtot),' relative: ',num2str(dstep/dtot)]);
    end
    
    %uses version of tidyconfig rewritten in C (tidyconfig scales by
    %n*#groups)
    [S,y] = tidyconfig_c(S,y);  %note tidyconfig reorders along node numbers
     
    
    %calculate modularity and return if converged
    if isequal(Sb,S)
        Q=0;
        P=sparse(y,1:length(y),1);
        for i=1:length(M(1))
            Q=Q+(P*M(i))'*P(:,i);
        end
        return
    end
        
    
    %check wether #groups < limit 
    t = length(unique(S));
    if (t>limit)
       M=@(i) metanetwork_i(B,S,t,i);   %use function handle if #groups>limit 
    else
        J = zeros(t);   %convert to matrix if #groups small enough
        for c=1:t
            J(:,c)=metanetwork_i(B,S,t,c);
        end
        B = J;
        M=B;
    end
    
end

    
    
    

S2 = (1:length(B))';
Sb = [];
while ~isequal(Sb,S2) %This is the loop around Blondel et al's "second phase"
    %keyboard
    y = unique(S2);  %unique also puts elements in ascending order
    Sb = S2;

    clocktime=clock;
    disp(['Merging ',num2str(length(y)),' communities  ',num2str(clocktime(4:6))]);
    yb = [];
    
    G = sparse(1:length(y),y,1);
    P=G';
    Q = sum(sum((P*M).*(P)));
    Qb = -inf;
    
    while (~isequal(yb,y)) && ((Q-Qb)/abs(Q)>2*eps) %This is the loop around Blondel et al's "first phase"
        
        disp([num2str(length(unique(y))),' ',num2str(Q)])
        yb = y;
        Qb=Q;
        
        for i = randperm(length(M))
            u = unique([y(i);y(M(:,i)>0)]);
           
            %dH = modchange_y(M(:,i),y,u); %relative changes in modularities
            
            dH = (M(:,i)'*G(:,u));
            
            dH(u==y(i)) = dH(u==y(i)) - M(i,i);
            [~, k] = max(dH);
            %%only move to different group if it is more optimized than
            %%staying in same group (up to error with double precision)
            if(dH(k)>(dH(u==y(i))+1*eps))
                y(i) = u(k);
                G(i,y(i))=0;
                G(i,u(k))=1;
            end
        end
        
        P=G';
        Q = sum(sum((P*M).*(P)));
        
    end
    
    y = tidyconfig(y);  %note tidyconfig reorders along node numbers
    for i = 1:length(y)
        S(S==i) = y(i);
        S2(S2==i) = y(i);
    end
   
    M = metanetwork(B,S2);    
end

%-----%
function M = metanetwork(J,S)
%Computes new aggregated network (communities --> nodes)
if(issparse(J))
    m=max(S);
    [i,j,v]=find(J);
    M = sparse(S(i),S(j),v,m,m);
else
    PP = sparse(1:length(S),S,1);
    M = PP'*J*PP;
end

function Mi = metanetwork_i(J,S,t,i) 
%ith column of metanetwork (used to create function handle)
%J is a function handle
  
    Mi=sparse([],[],[],t,1,t);
    for j=find(S==i)'
        Jj=J(j);
        [ii,k,v]=find(Jj);
        Mi=Mi+sparse(S(ii),k,v,t,1);       
    end
    
  
            

%-----%
function S = tidyconfig(S)
%This function remains almost identical to that originally written by
%Stephen Reid for his greedy.m code.
%   tidy up S i.e.  S = [2 4 2 6] -> S = [1 2 1 3]
T = zeros(length(S),1);
for i = 1:length(S)
    if T(i) == 0
        T(S==S(i)) = max(T) + 1;
    end
end
S = T;
