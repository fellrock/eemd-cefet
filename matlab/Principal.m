clc
clear all

M = 10;
g = modelos( 1 );
g = g( 1 : M );

N = 5;
P = 10;
numberOfIterations = 10000;
numberOfRepeats = 10;
beta = 1.5e-2;
betaG = 1e-3;
betaW = 10 * beta;
filtro = [ 1 -0.8 ];
sigmanu2 = 1e-4;
w0 = zeros( M + N - 1 , P );
sigmag2 = 1e-7;
actualG = g + sqrt( sigmag2 ) * randn( size( g ) );

for p = 1 : P
    
    w0( : , p ) = conv( actualG , randn( N , 1 ) );
    
end

% LMS
MSD{ 1 } = 0;
MSE{ 1 } = 0;

for p = 1 : P
   
    [ currentMSD , currentMSE ] = runLMS( w0( : , p ) , numberOfIterations , numberOfRepeats , beta , filtro , sigmanu2 );
    MSD{ 1 } = MSD{ 1 } + currentMSD / P;
    MSE{ 1 } = MSE{ 1 } + currentMSE / P;
    
end

[ MSD{ 2 } , MSE{ 2 } ] = runManifoldLMS( w0 , numberOfIterations , numberOfRepeats , betaG , betaW , filtro , sigmanu2 , N , M , g );


set( figure , 'Color' , 'w' )
plot( 1 : numberOfIterations , 10 * log10( MSD{ 1 } ) , 'b' )
hold on
plot( 1 : numberOfIterations , 10 * log10( MSD{ 2 } ) , 'r' )
xlabel( 'Iteration number' )
ylabel( 'MSD (dB)' )
grid on
axis tight
    
set( figure , 'Color' , 'w' )
plot( 1 : numberOfIterations , 10 * log10( MSE{ 1 } ) , 'b' )
hold on
plot( 1 : numberOfIterations , 10 * log10( MSE{ 2 } ) , 'r' )
xlabel( 'Iteration number' )
ylabel( 'MSE (dB)' )
grid on
axis tight
