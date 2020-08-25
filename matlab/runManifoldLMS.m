function [ MSD , MSE ] = runManifoldLMS( w0 , numberOfIterations , numberOfRepeats , betaG , betaW , filtro , sigmanu2 , N , M , g0 )

MSD = zeros( numberOfIterations , 1 );
MSE = zeros( numberOfIterations , 1 );
P = size( w0 , 2 );

x  = cell( P , 1 );
d  = cell( P , 1 );
xk = cell( P , 1 );
f  = zeros( M + N - 1 , P );
ek = zeros( P , 1 );
Wk = cell( P );

for p = 1 : P
    
    Wk{ p } = zeros( M + N - 1 , M );
    
end

for repeat = 1 : numberOfRepeats
    
    for p = 1 : P
        
        x{ p } = randn( numberOfIterations + M + N - 1 , 1 );
        x{ p } = filter( filtro , 1 , x{ p } );
        d{ p } = filter( w0( : , p ) , 1 , x{ p } );
        d{ p } = d{ p } + sqrt( sigmanu2 ) * randn( size( d{ p } ) );
        
    end
    
    gk = g0;%zeros( M , 1 ) + 1e-12;
    wk = zeros( N , P ) + 1e-12;
    
    for k = N + M - 1 : numberOfIterations + N + M - 2
        
        Deltag = 0;
        for p = 1 : P
        
            xk{ p } = x{ p }( k : -1 : k - M - N + 2 );
            f( : , p ) = conv( gk , wk( : , p ) );
            ek( p ) = d{ p }( k ) - transpose( f( : , p ) ) * xk{ p };
            
            Wk{ p }( : , : ) = 0;
            wE = zeros( 2 * N - 2 + M , 1 );
            wE( 1 : N ) = wk( end : -1 : 1 , p ); 
            
            for n = 1 : M + N - 1
                
                Wk{ p }( n , : ) = wE( N : N + M - 1 );
                wE( 2 : end ) = wE( 1 : end - 1 );
                wE( 1 ) = 0;
                
            end
            
            Deltag = Deltag + betaG * ek( p ) * transpose( Wk{ p } ) * xk{ p };
            
            for n = 0 : N - 1
                
                wk( n + 1 , p ) = wk( n + 1 , p ) + betaW * ek( p ) * transpose( x{ p }( k - n : -1 : k - n - M + 1 ) ) * gk; 
                
            end
            
        end
        
        gk = gk + Deltag;
%        gk = gk / norm( gk );
        
        for p = 1 : P
            
            f( : , p ) = conv( gk , wk( : , p ) );
            MSD( k - N - M + 2 ) = MSD( k - N - M + 2 ) + norm( f( : , p ) - w0( : , p ) ) ^ 2 / numberOfRepeats / P;  
    
        end
                
        MSE( k - N - M + 2 ) = MSE( k - N - M + 2 ) + mean( ek .^ 2 ) / numberOfRepeats ;  
        
    end
    
end

