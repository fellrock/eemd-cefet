function [ MSD , MSE ] = runLMS( w0 , numberOfIterations , numberOfRepeats , beta , filtro , sigmanu2 )

N = numel( w0 );

MSD = zeros( numberOfIterations , 1 );
MSE = zeros( numberOfIterations , 1 );

for repeat = 1 : numberOfRepeats
    
    x = randn( numberOfIterations + N - 1 , 1 );
    x = filter( filtro , 1 , x );
    d = filter( w0 , 1 , x );
    d = d + sqrt( sigmanu2 ) * randn( size( d ) );
    wk = zeros( N , 1 );
    
    for k = N : numberOfIterations + N - 1
        
        xk = x( k : -1 : k - N + 1 );
        yk = transpose( wk ) * xk;
        ek = d( k ) - yk;
        wk = wk + beta * ek * xk;
        
        MSD( k - N + 1 ) = MSD( k - N + 1 ) + norm( wk - w0 ) ^ 2 / numberOfRepeats;  
        MSE( k - N + 1 ) = MSE( k - N + 1 ) + ek ^ 2 / numberOfRepeats;  
        
    end
    
end

