function DodecahedronSimulation
% load('errTotal1.000e-05.mat')
% errorVecMat
% return

n=1;

% % load('HybrideFullAllowedLookupZerr');
% % 
% % newMatZHybride
% 
% SynZ1 = zeros(1,30);
% SynZ1(1,1) = 1;
% SynZ1(1,4) = 1;
% SynZ2 = zeros(1,30);
% SynZ2(1,1) = 1;
% SynZ2(1,4) = 1;
% SynZ3 = zeros(1,30);
% SynZ3(1,1) = 1;
% SynZ3(1,4) = 1;
% 
% 
% tic
% CorrectionX(SynZ1,SynZ2,SynZ3)
% toc
% 
% return


numIterations = 10^7;
for i = 1:12
    v = [9*10^-6,10^-5,2*10^-5,3*10^-5,4*10^-5,5*10^-5,6*10^-5,7*10^-5,8*10^-5,9*10^-5,10^-4,2*10^-4];
    
    errRate = v(1,i);
    str_errRate = num2str(errRate,'%0.3e.mat');
    errX = 'errX';
    errZ = 'errZ';
    errY = 'errY';
    errTotal = 'errTotal';
    str_errX = strcat(errX,str_errRate);
    str_errY = strcat(errY,str_errRate);
    str_errZ = strcat(errZ,str_errRate);
    str_errTotal = strcat(errTotal,str_errRate);
    [logX,logZ,logY,logTotal] = depolarizingSimulator(numIterations,errRate,n);
    parsaveErrorVec(str_errX,logX);
    parsaveErrorVec(str_errY,logY);
    parsaveErrorVec(str_errZ,logZ);
    parsaveErrorVec(str_errTotal,logTotal);
end


end

function parsaveErrorVec(fname,errorVecMat)
save(fname,'errorVecMat');
end

function Output = PropagationStatePrepArb(C, n, e)

% C encodes information about the circuit
% n is the number of qubits per encoded codeblock
% eN are the error entries (location and time)

% The Output matrix will contain (2q+m) rows and n columns. The paramters Output(i,j) are:
% i <= 2q and odd: Stores X type errors for output qubit #ceil(i/2) (order based on matrix C)
% i <= 2q and even: Stores Z type errors for output qubit #ceil(i/2) (order based on matrix C)
% i > 2q: Stores measurement error for measurement number #i-2q (type is X or Z depending on measurement type)

% The Errors matrix will be a tensor with 2 parameters Output(i, j)
% i: logical qubit number
% j<=n: X error at physical qubit number j
% j>n : Z error at physical qubit number j


% Error inputs eN are characterized by four parameters, thus a vector (i,j, k, l)
% i: error type: 0 = I,  1 = X, 2 = Z, 3 = Y
% j: logical qubit number (if the entry here is '0' then this will correspond to an identity "error")
% k: physical qubit number within a codeblock
% l: time location

N_meas = length(find(C==5)) + length(find(C==6)); % number of logical measurements in the circuit
N_output = length(C(:,end)) - sum(C(:,end)==-1) - sum(C(:,end)==5) - sum(C(:,end)==6); % number of output active qubits that have not been measured

Meas_dict = [find(C==5); find(C==6)];
Meas_dict = sort(Meas_dict);

Output = zeros(2*N_output + N_meas, n);
Errors = zeros(length(C(:,1)), 2*n);


for t= 1:length(C(1,:))
    
    % If the error occurs at a measurement location then the error is introduced before propagation of faults
    for j = 1:length(e(:,1))
        if e(j,1) ~= 0
            if e(j,4) == t && ( ( C(e(j,2),t) == 5 ) || ( C(e(j,2),t) == 6 ) )
                if e(j,1) == 1
                    Errors(e(j,2), e(j,3)) = mod(Errors(e(j,2), e(j,3)) + 1, 2);
                end
                if e(j,1) == 2
                    Errors(e(j,2), e(j,3)+n) = mod(Errors(e(j,2), e(j,3)+n) + 1, 2);
                end
                if e(j,1) == 3
                    Errors(e(j,2), e(j,3)) = mod(Errors(e(j,2), e(j,3)) + 1, 2);
                    Errors(e(j,2), e(j,3)+n) = mod(Errors(e(j,2), e(j,3)+n) + 1, 2);
                end
            end
        end
    end
    
    % Propagation of errors (do not need to do it for the first step of circuit
    if t>1
        for i = 1:length(C(:,t))
            if C(i, t) == 10
                % In this case must flip the X and Z error information
                v1 = Errors(i, 1:n);
                v2 = Errors(i, (n+1):end);
                Errors(i, 1:n) = v2;
                Errors(i, n+1:end) = v1;
            end
            
            if C(i, t) == 11
                v1 = Errors(i, 1:n);
                v2 = Errors(i, (n+1):end);
                Errors(i, 1:n) = mod(v1+v2, 2);
            end
            
            if C(i,t) > 1000
                % check to see if target qubit according to control qubit is actually a target qubit
                if C(C(i,t) - 1000, t) == 1000
                    v1 = Errors(i, 1:n);
                    v2 = Errors(i, (n+1):end);
                    w1 = Errors(C(i,t) - 1000, 1:n);
                    w2 = Errors(C(i,t) - 1000, (n+1):end);
                    %mod(v2+w2,2)
                    Errors(C(i,t) - 1000, 1:n) = mod(v1+w1, 2);
                    Errors(i, (n+1):end) = mod(v2+w2, 2);
                end
            end
            
            if C(i,t) == 5
                % This corresponds to X measurement, therefore need to look at Z errors
                find(Meas_dict==(i+(t-1)*length(C(:,1))) ); % Dont think these are needed, used for checking
                Output(2*N_output + find(Meas_dict==(i+(t-1)*length(C(:,1))) ), :) = Errors(i, (n+1):end);
                Errors(i, 1:end) = 0;
            end
            
            if C(i,t) == 6
                % This corresponds to Z measurement, therefore need to look at Z errors
                find(Meas_dict==(i+(t-1)*length(C(:,1))) ); % Dont think these are needed, used for checking
                Output(2*N_output + find(Meas_dict==(i+(t-1)*length(C(:,1))) ), :) = Errors(i, 1:n);
                Errors(i, 1:end) = 0;
            end
            
        end
    end
    
    % Introduce faults for locations that are not measurements
    for j = 1:length(e(:,1))
        if e(j,1) ~= 0
            if e(j,4) == t && ( C(e(j,2),t) ~= 5 ) && ( C(e(j,2),t) ~= 6 )
                % This IF statement checks to see if the gate at this location is NOT a CNOT or Prep
                if ( C(e(j,2),t) < 1000 ) %&& ( C(e(j,2),t) ~= 3 ) && ( C(e(j,2),t) ~= 4 )
                    if e(j,1) == 1
                        Errors(e(j,2), e(j,3)) = mod(Errors(e(j,2), e(j,3)) + 1, 2);
                    end
                    if e(j,1) == 2
                        Errors(e(j,2), e(j,3)+n) = mod(Errors(e(j,2), e(j,3)+n) + 1, 2);
                    end
                    if e(j,1) == 3
                        Errors(e(j,2), e(j,3)) = mod(Errors(e(j,2), e(j,3)) + 1, 2);
                        Errors(e(j,2), e(j,3)+n) = mod(Errors(e(j,2), e(j,3)+n) + 1, 2);
                    end
                end
                % Introduce errors in the case of CNOT gate for control and target qubits
                % Errors for control qubit are entry mod(e(j,1),4) according to standard indexing above
                if ( C(e(j,2),t) > 1000 )
                    if C(C(e(j,2),t) - 1000, t) == 1000
                        if mod(e(j,1),2) == 1
                            Errors(e(j,2), e(j,3)) = mod(Errors(e(j,2), e(j,3)) + 1, 2);
                        end
                        if mod(e(j,1),4) > 1
                            Errors(e(j,2), e(j,3)+n) = mod(Errors(e(j,2), e(j,3)+n) + 1, 2);
                        end
                        if mod(floor(e(j,1)/4),2) == 1
                            Errors(C(e(j,2),t) - 1000, e(j,3)) = mod(Errors(C(e(j,2),t) - 1000, e(j,3)) + 1, 2);
                        end
                        if mod(floor(e(j,1)/4),4) > 1
                            Errors(C(e(j,2),t) - 1000, e(j,3)+n) = mod(Errors(C(e(j,2),t) - 1000, e(j,3)+n) + 1, 2);
                        end
                    end
                end
            end
        end
    end
    
end

%Errors
counter = 1; % This will be used to iterate over the different qubits
for j = 1:length(C(:,end))
    if (C(j,end) ~= -1) && (C(j,end) ~= 5) &&  (C(j,end) ~= 6)
        Output(counter,:) = Errors(j,1:n);
        Output(counter+1,:) = Errors(j,(n+1):end);
        counter = counter + 2;
    end
end

end

function Output = ErrorGenerator(Cmat,errRate)
%This function outputs an error vector based on the input circuit
%represented by Cmat. We have the following noise model.

% 1) |0> state preparation: Perfect |0> state followed by X error with probability p
% 2) Z-measurement: X pauli with probability p followed by perfect Z-basis
% measurement.
% 3) CNOT: Perfect CNOT followed by
%{IX,IY,IZ,XI,YI,ZI,XX,XY,XZ,ZX,ZY,ZZ,YX,YY,YZ} with probability p/15 each.
% 4) Hadamard: Perfect Hadamard followed by {X,Y,Z} with probability p/12
%each.
% 5) SWAP: Perfect SWAP followed by
% {IX,IY,IZ,XI,YI,ZI,XX,XY,XZ,ZX,ZY,ZZ,YX,YY,YZ} with probability p/60
% each.
% Storage: Pauli {X,Y,Z} with probability p/30 each.

% Here errRate = p and Cmat is the circuit representing the surface code
% lattice.


e = zeros(1,4);
counter = 1;

for i = 1:length(Cmat(:,1))
    for j = 1:length(Cmat(1,:))
        
        % Adds storage errors with probability p/10
        if (Cmat(i,j) == 1)
            xi = rand;
            if xi < errRate/10
                k = randi([1,3]);
                e(counter,:) = [k,i,1,j];
                counter = counter + 1;
            end
        end
        
        % Adds state-preparation errors with probability 2p/3
        if Cmat(i,j) == 4
            xi = rand;
            if xi < 2*errRate/3
                e(counter,:) = [1,i,1,j];
                counter = counter + 1;
            end
        end
        if Cmat(i,j) == 3
            xi = rand;
            if xi < 2*errRate/3
                e(counter,:) = [2,i,1,j];
                counter = counter + 1;
            end
        end
        
        % Adds measurement errors with probability 2p/3
        if Cmat(i,j) == 6
            xi = rand;
            if xi < 2*errRate/3
                e(counter,:) = [1,i,1,j];
                counter = counter + 1;
            end
        end
        if Cmat(i,j) == 5
            xi = rand;
            if xi < 2*errRate/3
                e(counter,:) = [2,i,1,j];
                counter = counter + 1;
            end
        end
        
        
        % Adds CNOT errors with probability p
        if (Cmat(i,j) > 1000)
            xi = rand;
            if xi < errRate
                k = randi([1,15]);
                e(counter,:) = [k,i,1,j];
                counter = counter + 1;
            end
        end
        
    end
end

Output = e;

end

function Output = CorrectionX(SynX1,SynX2,SynX3)

% This function uses our decoding algorithm to correct X errors based on
% the measured syndromes SynX1, SynX2 and SynX3 during rounds 1, 2 and 3.

HybrideLookupTableX = [0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     1     0     0     0
     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     1     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0];

s1 = max(SynX1);
s2 = max(SynX2);
s3 = max(SynX3);

if (s1 == 0 && (s2 == 0 || s3 == 0)) || (s2 == 0 && s3 == 0) % If at least two syndromes are trivial, apply no correction
    corrX = zeros(1,30);
elseif isequal(SynX1,SynX2) || isequal(SynX1,SynX3)  % If at least two syndromes are equal
    correctionRow = 1;
    for ii = 1:length(SynX1)
        correctionRow = correctionRow + 2^(11-ii)*SynX1(ii);
    end
    corrX = HybrideLookupTableX(correctionRow,:);
elseif isequal(SynX2,SynX3)
    correctionRow = 1;
    for ii = 1:length(SynX2)
        correctionRow = correctionRow + 2^(11-ii)*SynX2(ii);
    end
    corrX = HybrideLookupTableX(correctionRow,:);
else % If all three syndromes are non-trivial and different, use the last syndrome to correct
    correctionRow = 1;
    for ii = 1:length(SynX3)
        correctionRow = correctionRow + 2^(11-ii)*SynX3(ii);
    end
    corrX = HybrideLookupTableX(correctionRow,:);
end

Output = corrX;

end

function Output = CorrectionZ(SynZ1,SynZ2,SynZ3)

% This function uses our decoding algorithm to correct Z errors based on
% the measured syndromes SynZ1, SynZ2 and SynZ3 during rounds 1, 2 and 3.

HybrideLookupTableZ = [0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     1     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0
     1     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     1     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1
     1     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1     0     1     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     1     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     1
     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     1     0
     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     1     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     1     0     0
     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     0     1     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
     0     0     0     0     1     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     0     1     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     0     1     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     0     1     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1
     1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     1     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     1
     1     0     0     0     0     0     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     1     0
     1     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     1     0
     1     0     0     0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     1     0     1];

s1 = max(SynZ1);
s2 = max(SynZ2);
s3 = max(SynZ3);

if (s1 == 0 && (s2 == 0 || s3 == 0)) || (s2 == 0 && s3 == 0) % If at least two syndromes are trivial, apply no correction
    corrZ = zeros(1,30);
elseif isequal(SynZ1,SynZ2) || isequal(SynZ1,SynZ3)  % If at least two syndromes are equal
    correctionRow = 1;
    for ii = 1:length(SynZ1)
        correctionRow = correctionRow + 2^(11-ii)*SynZ1(ii);
    end
    corrZ = HybrideLookupTableZ(correctionRow,:);
elseif isequal(SynZ2,SynZ3)
    correctionRow = 1;
    for ii = 1:length(SynZ2)
        correctionRow = correctionRow + 2^(11-ii)*SynZ2(ii);
    end
    corrZ = HybrideLookupTableZ(correctionRow,:);
else % If all three syndromes are non-trivial and different, use the last syndrome to correct
    correctionRow = 1;
    for ii = 1:length(SynZ3)
        correctionRow = correctionRow + 2^(11-ii)*SynZ3(ii);
    end
    corrZ = HybrideLookupTableZ(correctionRow,:);
end

Output = corrZ;

end

function Output = ConvertErrorXStringToErrorVec(ErrStrX)

e = zeros(1,4);
counter = 1;

for i = 1:length(ErrStrX(1,:))
    if ErrStrX(1,i) == 1
        e(counter,:) = [1,i+11,1,1];
        counter = counter + 1;
    end
end

Output = e;

end

function Output = ConvertErrorZStringToErrorVec(ErrStrZ)

e = zeros(1,4);
counter = 1;

for i = 1:length(ErrStrZ(1,:))
    if ErrStrZ(1,i) == 1
        e(counter,:) = [2,i+11,1,1];
        counter = counter + 1;
    end
end

Output = e;

end

function [OutputX,OutputZ,OutputY,OutputTotal] = depolarizingSimulator(numIterations,errRate,n)

% This function performs numIterations iterations of the depolarizing
% channel applied to the [[30,8,3]] dodecahedron code. 

logicalX = 0;
logicalZ = 0;
logicalY = 0;
logicalTotal = 0;

% Circuit descriptors:
% -1: qubit non-active
% 0: Noiseless memory
% 1: Gate memory
% 2: Measurement memory
% 3: Preparation in X basis (|+> state)
% 4: Preparation in Z basis (|0> state)
% 5: Measurement in X basis
% 6: Measurement in Z basis
% 7: X gate
% 8: Z gate
% 10: H gate
% 11: S gate
% 20: T gate
%1---: Control qubit for CNOT with target ---
%1000: Target qubit for CNOT

load('HybrideFullAllowedLookupXerr');
load('HybrideFullAllowedLookupZerr');

% X stabilizer generators
gX = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0
    1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0
    0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0
    0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1
    0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0
    0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1];

% X logical operators
XL = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
    1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1];

% Z stabilizer generators
gZ = [0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
    0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0
    1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0
    0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,1,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,0,0
    0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0
    0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0];

% Z logical operators
ZL = [1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0
    0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
    0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1
    0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0
    0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0
    0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0];


% Circuit representing one round of stabilizer measurement
Circuit = [3, 1013, 1012, 1014, 1016, 1015, 0, 0, 0, 0, 0, 5
    3, 1017, 1018, 1020, 1019, 1021, 0, 0, 0, 0, 0, 5
    3, 1022, 1025, 1024, 1023, 1026, 0, 0, 0, 0, 0, 5
    3, 1027, 1017, 1028, 1030, 1029, 0, 0, 0, 0, 0, 5
    3, 1018, 1032, 1022, 1031, 1033, 0, 0, 0, 0, 0, 5
    3, 1034, 1023, 1027, 1035, 1036, 0, 0, 0, 0, 0, 5
    3, 1012, 1024, 1037, 1034, 1038, 0, 0, 0, 0, 0, 5
    3, 1019, 1013, 1040, 1039, 1028, 0, 0, 0, 0, 0, 5
    3, 1041, 1014, 1031, 1025, 1037, 0, 0, 0, 0, 0, 5
    3, 1035, 1015, 1038, 1029, 1039, 0, 0, 0, 0, 0, 5
    3, 1040, 1020, 1016, 1032, 1041, 0, 0, 0, 0, 0, 5
    1, 1000, 1000, 1, 1, 1, 1046, 1, 1, 1, 1045, 1
    1, 1000, 1000, 1, 1, 1, 1, 1047, 1, 1, 1046, 1
    1, 1, 1000, 1000, 1, 1, 1, 1043, 1, 1, 1047, 1
    1, 1, 1000, 1, 1, 1000, 1, 1044, 1, 1043, 1, 1
    1, 1, 1, 1000, 1000, 1, 1, 1, 1045, 1044, 1, 1
    1, 1000, 1000, 1, 1, 1, 1044, 1042, 1, 1, 1, 1
    1, 1000, 1000, 1, 1, 1, 1, 1, 1047, 1042, 1, 1
    1, 1000, 1, 1, 1000, 1, 1047, 1, 1050, 1, 1, 1
    1, 1, 1000, 1000, 1, 1, 1051, 1, 1044, 1, 1, 1
    1, 1, 1, 1, 1, 1000, 1, 1051, 1, 1050, 1, 1
    1, 1000, 1, 1000, 1, 1, 1, 1045, 1, 1, 1042, 1
    1, 1, 1000, 1, 1000, 1, 1042, 1, 1, 1, 1043, 1
    1, 1, 1000, 1000, 1, 1, 1045, 1, 1, 1052, 1, 1
    1, 1, 1000, 1, 1000, 1, 1, 1, 1043, 1051, 1, 1
    1, 1, 1, 1, 1, 1000, 1052, 1, 1051, 1, 1, 1
    1, 1000, 1, 1000, 1, 1, 1, 1046, 1042, 1, 1, 1
    1, 1, 1, 1000, 1, 1000, 1048, 1, 1, 1046, 1, 1
    1, 1, 1, 1, 1000, 1000, 1, 1052, 1, 1, 1044, 1
    1, 1, 1, 1, 1000, 1, 1, 1048, 1052, 1, 1, 1
    1, 1, 1, 1000, 1000, 1, 1, 1049, 1, 1047, 1, 1
    1, 1, 1000, 1, 1000, 1, 1, 1, 1048, 1045, 1, 1
    1, 1, 1, 1, 1, 1000, 1, 1, 1, 1049, 1048, 1
    1, 1000, 1, 1, 1000, 1, 1049, 1, 1046, 1, 1, 1
    1, 1000, 1, 1, 1000, 1, 1043, 1050, 1, 1, 1, 1
    1, 1, 1, 1, 1, 1000, 1050, 1, 1049, 1, 1, 1
    1, 1, 1, 1000, 1, 1000, 1, 1, 1, 1, 1049, 1
    1, 1, 1, 1000, 1, 1000, 1, 1, 1, 1, 1052, 1
    1, 1, 1, 1, 1000, 1000, 1, 1, 1, 1, 1050, 1
    1, 1000, 1, 1000, 1, 1, 1, 1, 1, 1048, 1, 1
    1, 1000, 1, 1, 1, 1000, 1, 1, 1, 1, 1051, 1
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6
    4, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000, 6];

parfor iii = 1:numIterations
    e = zeros(1,4);
    XSyn = zeros(3,11);
    ZSyn = zeros(3,11);
    for i = 1:3
        eCircuit = ErrorGenerator(Circuit,errRate);
        eTemp = [e;eCircuit];
        Cout = transpose(PropagationStatePrepArb(Circuit, n, eTemp));
        Xerror = Cout(1,1:2:60);
        Zerror = Cout(1,2:2:60);
        eX = ConvertErrorXStringToErrorVec(Xerror);
        eZ = ConvertErrorZStringToErrorVec(Zerror);
        ZSyn(i,:) = Cout(1,61:71);
        XSyn(i,:) = Cout(1,72:82);
        if i ~= 3
            e = [eX;eZ];
        end
    end
    
    SynX1 = XSyn(1,:);
    SynX2 = XSyn(2,:);
    SynX3 = XSyn(3,:);
    CorrX = CorrectionX(SynX1,SynX2,SynX3);    
    eXfinal = mod(Xerror+CorrX,2);
    
    SynZ1 = ZSyn(1,:);
    SynZ2 = ZSyn(2,:);
    SynZ3 = ZSyn(3,:);
    CorrZ = CorrectionZ(SynZ1,SynZ2,SynZ3);
    eZfinal = mod(Zerror+CorrZ,2);
            
    eX = ConvertErrorXStringToErrorVec(eXfinal);        
    eZ = ConvertErrorZStringToErrorVec(eZfinal);    
    e = [eX;eZ];
    
    XSyn = zeros(3,11);
    ZSyn = zeros(3,11);
    for i = 1:3
        eCircuit = ErrorGenerator(Circuit,errRate);           
        eTemp = [e;eCircuit];
        Cout = transpose(PropagationStatePrepArb(Circuit, n, eTemp));
        Xerror = Cout(1,1:2:60);
        Zerror = Cout(1,2:2:60);
        eX = ConvertErrorXStringToErrorVec(Xerror);        
        eZ = ConvertErrorZStringToErrorVec(Zerror);        
        ZSyn(i,:) = Cout(1,61:71);
        XSyn(i,:) = Cout(1,72:82);
        e = [eX;eZ];
    end
    
    SynX1 = XSyn(1,:);
    SynX2 = XSyn(2,:);
    SynX3 = XSyn(3,:);
    CorrX = CorrectionX(SynX1,SynX2,SynX3);
    
    eXfinal = mod(Xerror+CorrX,2);
    
    SynZ1 = ZSyn(1,:);
    SynZ2 = ZSyn(2,:);
    SynZ3 = ZSyn(3,:);
    CorrZ = CorrectionZ(SynZ1,SynZ2,SynZ3);
    
    eZfinal = mod(Zerror+CorrZ,2);
    
    % Next we do one round of perfect error correction in order to determine if
    % there is a logical X or Z fault
    
    Syn = zeros(1,11);
    for i = 1:length(gZ(:,1))
        Syn(1,i) = mod(sum(conj(eXfinal).* gZ(i,:)),2);
    end
    
    correctionRow = 1;
    for i = 1:length(Syn)
        correctionRow = correctionRow + 2^(11-i)*Syn(i);
    end
    eXperfect = mod(newMatXHybride(correctionRow,:)+eXfinal,2);
    
    
    Syn = zeros(1,11);
    for i = 1:length(gX(:,1))
        Syn(1,i) = mod(sum(conj(eZfinal).* gX(i,:)),2);
    end
    
    correctionRow = 1;
    for i = 1:length(Syn)
        correctionRow = correctionRow + 2^(11-i)*Syn(i);
    end
    eZperfect = mod(newMatZHybride(correctionRow,:)+eZfinal,2);
    
    % Check if there is a logical fault by looking if eXperfect anti-commutes
    % with a Z logical operator and eZperfect anti-commuteswith an X logical
    % operator
    
    SynLogicalX = zeros(1,length(ZL(:,1)));
    for i = 1:length(ZL(:,1))
        SynLogicalX(1,i) = mod(sum(conj(eXperfect).* ZL(i,:)),2);
    end
    
    SynLogicalZ = zeros(1,length(XL(:,1)));
    for i = 1:length(XL(:,1))
        SynLogicalZ(1,i) = mod(sum(conj(eZperfect).* XL(i,:)),2);
    end
    
    
    if max(SynLogicalX) ~= 0
        logicalX = logicalX + 1;        
    end
    if max(SynLogicalZ) ~= 0
        logicalZ = logicalZ + 1;        
    end
    if (max(SynLogicalZ) ~= 0) && (max(SynLogicalX) ~= 0)
        logicalY = logicalY + 1;        
    end
    
    if max(SynLogicalX) ~= 0 || max(SynLogicalZ) ~= 0
        logicalTotal = logicalTotal + 1;
    end
    
end



OutputX = logicalX/numIterations;
OutputZ = logicalZ/numIterations;
OutputY = logicalY/numIterations;
OutputTotal = logicalTotal/numIterations;

end