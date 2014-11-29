
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

// DEFINIZIONE DI COSTANTI
#define N 100
#define ETA 0.05
#define PI 3.14159265358979

// PROTOTIPI DELLE FUNZIONI
double speed_gen( void );
double coll_time_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s );
double v_min_search( double v[ 9 ]);
double amongst_copies_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s );
void evolve( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], double delta );
void speed_refresh( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], int p1, int p2 ); // p1 e p2 sono risp. le posizioni delle due particelle che collidono

int main( void ){
	// DICHIARAZIONE DELLE VARIABILI

	int n = sqrt( N );

	// variabili per l'inizializzazione
	double	sigma = pow( ( 4.0 * ETA )/( PI * N ) , 0.5 ); // diametro dei dischi
	double v0x[N], v0y[N], r0x[N], r0y[N]; // posizioni e velocità iniziali
	double Mtime[N][N]; // matrice dei tempi
	double step = 1/(double)n; // distanza fra i centri di dischi nella configurazione iniziale
	double vcmx=0, vcmy=0; // velocità del centro di massa
	double ctime, min;

	// variabili per l'evoluzione
	double t0 = 0; // tempo iniziale al passo i-esimo
	double delta_t; // intervallo fra due passi 
	double vx[ N ], vy[ N ], rx[ N ], ry[ N ];

	// indici per cicli
	int i, j, l, m, min_place[ 2 ];

	// gestione file
	FILE *cfPtr;
	cfPtr = fopen( "data.dat", "w" );

	srand( time( NULL ) );// seme della funzione random
	
	/* printf( "\n#%g#\n", sigma ); */

	// coordinate x dei punti
	for( i=0; i<N; i++ ){
		r0x[i] = ( ( i % n ) * step );
	}	
	// coordinate y dei punti
	for( i=0; i<N; i++ ){
		r0y[i] = ( i / n ) * step;
	}

	// assegno le velocità a caso
	for( i=0; i<N; i++ ){
		v0x[i] = speed_gen();
		v0y[i] = speed_gen();
	}
	
	// calcolo la velocità del centro di massa
	for( i=0; i<N; i++ ){
		vcmx += (double)v0x[i] / N;
		vcmy += (double)v0y[i] / N;
	}
	// riassegno le velocità con la condizione vcm=0
	for( i=0; i<N; i++ ){
		v0x[i] = v0x[i] - vcmx;
		v0y[i] = v0y[i] - vcmy;
	}
	
	vcmx=0, vcmy=0;
	// controllo che il centro di massa sia fermo
		for( i=0; i<N; i++ ){
		vcmx += (double)v0x[i] / N;
		vcmy += (double)v0y[i] / N;
	}
	
	/*
	printf( "\n Esempio= %f\n" , amongst_copies_search( r0x[1], r0y[1], v0x[1], v0y[1], r0x[0], r0y[0], v0x[0], v0y[0], sigma ) );
	// controlli vari
	printf( "%f,%f\n\n", vcmx, vcmy );
	for( i=0; i<N; i++ ){
		printf( "posizioni:\t(%f,%f)", r0x[i], r0y[i] );
		printf( "\t velocità:\t(%f,%f)", v0x[i], v0y[i] );
		printf( "\n\n" );
	}
	*/
		
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			ctime = amongst_copies_search( r0x[i], r0y[i], v0x[i], v0y[j], r0x[j], r0y[j], v0x[j], v0y[j], sigma );
			if( i==j ) Mtime[i][j] = 0;
			else if( ctime == FLT_MAX ) Mtime[i][j] = INFINITY;
			else Mtime[i][j] = ctime;
		}
	}
	
	/*
	printf( "\n Tabella\n -------------------------------------------\n" );
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			printf( "\t%.2f", Mtime[i][j] );
		}
		printf( "\n" );
	}
	printf( "\n -------------------------------------------" );
	*/
	
	// determino il tempo di collisione più piccolo


	/***************************** algoritmo di scansione completa della matrice dei tempi
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			if( Mtime[i][j] < min && i != j ){
				min = Mtime[ i ][ j ];
				min_place[ 0 ] = i;
				min_place[ 1 ] = j;
			}
		}
	}
	printf( "\n min = %g.\n", min );
	*/


	min = INFINITY;
	for( i=0; i<N; i++ ){
		for( j=i; j<N; j++ ){
			if( Mtime[i][j] < min && i != j ){
				min = Mtime[ i ][ j ];
				min_place[ 0 ] = i;	
				min_place[ 1 ] = j;
			}
		}
	}

	printf( "Il disco %d collide col disco %d al tempo t=%E.\n", min_place[ 0 ], min_place[ 1 ], min );
	
	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t x1\t y1\t x2\t y2 \t vx1\t vy1\t vx2\t vy2\n" );
		fprintf( cfPtr, "\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n", r0x[ min_place[ 0 ] ], r0y[ min_place[ 0 ] ], r0x[ min_place[ 1 ] ], r0y[ min_place[ 1 ] ], v0x[ min_place[ 0 ] ], v0y[ min_place[ 0 ] ], v0x[ min_place[ 1 ] ], v0y[ min_place[ 1 ] ] );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono:\t %f.\n",  pow( r0x[ min_place[ 0 ] ] - r0x[ min_place[ 1 ] ], 2 ) + pow(r0y[ min_place [ 0 ] ] - r0y[ min_place[ 1 ] ], 2 ) - sigma );
		} // chiudo else - scrittura file
	fclose( cfPtr );
		
	delta_t = min;
	for( i=0; i<N; i++ ){
		rx[ i ] = r0x[ i ];
		ry[ i ] = r0y[ i ];
		vx[ i ] = v0x[ i ];
		vy[ i ] = v0y[ i ];
	}

	// determino la posizione della particella 2 coinvolta nell'urto
	if( abs( rx[ min_place[ 0 ] ] - rx[ min_place[ 1 ] ]) > 0.5 )
		rx[ min_place[ 1  ] ] += floor( 0.5 +  rx[ min_place[ 0 ] ] - rx[ min_place[ 1 ] ] );

 
	if( abs( ry[ min_place[ 0 ] ] - ry[ min_place[ 1 ] ]) > 0.5 )
		ry[ min_place[ 1  ] ] += floor( 0.5 +  ry[ min_place[ 0 ] ] - ry[ min_place[ 1 ] ] );
	
	cfPtr = fopen( "data.dat", "a" );
	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono dopo l'aggiustamento:\t %f.\n",  pow( rx[ min_place[ 0 ] ] - rx[ min_place[ 1 ] ], 2 ) + pow(ry[ min_place [ 0 ] ] - ry[ min_place[ 1 ] ], 2 ) - sigma );

	// evolve( rx, ry, vx, vy, delta_t );
	
	// aggiorno le velocità delle paricelle che collidono
	speed_refresh( rx, ry, vx, vy,  min_place[ 0 ], min_place[ 1 ] );	

	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t x1\t y1\t x2\t y2 \t vx1\t vy1\t vx2\t vy2\n" );
		fprintf( cfPtr, "\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n", rx[ min_place[ 0 ] ], ry[ min_place[ 0 ] ], rx[ min_place[ 1 ] ], ry[ min_place[ 1 ] ], vx[ min_place[ 0 ] ], vy[ min_place[ 0 ] ], vx[ min_place[ 1 ] ], vy[ min_place[ 1 ] ] );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono:\t %f.\n",  pow( rx[ min_place[ 0 ] ] - rx[ min_place[ 1 ] ], 2 ) + pow(ry[ min_place [ 0 ] ] - ry[ min_place[ 1 ] ], 2 ) - sigma );
		} // chiudo else - scrittura file
	fclose( cfPtr );
		

	return 0;
}

// DICHIARAZIONE DELLE FUNZIONI

// genera le velocità casualmente
double speed_gen( void ){
	return ( rand()/(double)RAND_MAX ) * pow( ( -1 ), ( rand() % 2 ) );
}

// determina il tempo di sollisione fra due dischi
double coll_time_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s ){
	double dr_x, dr_y, dv_x, dv_y, dot_prod, sq_norm_r, sq_norm_v;
	
	// calcoli intermedi
	dr_x = r1x - r2x;
	dr_y = r1y - r2y;
	dv_x = v1x - v2x;
	dv_y = v1y - v2y;
	
	dot_prod = dr_x * dv_x + dr_y * dv_y;
	sq_norm_r = dr_x * dr_x + dr_y * dr_y;
	sq_norm_v = dv_x * dv_x + dv_y * dv_y;
	
	// uscita della funzione
	/*
	printf( "\n Cerca: %g, %g, %g, %e, %g\n", dot_prod, sq_norm_r, sq_norm_v, s * s, sq_norm_r-s*s );
	printf( " tempo: \%g\n", (- dot_prod - sqrt( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) ) ) / sq_norm_v );
	*/
	
	if( dot_prod < 0 && ( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) >=0 ) ) // controllo che i dischi si avvicinino e che i tempi siano reali
		return ( - dot_prod - sqrt( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) ) ) / sq_norm_v;
	else return FLT_MAX;
	
}

// determina il minimo di un vettore
double v_min_search( double v[ 9 ]){
	int k=0;
	double min = v[0];
	
	for( k=1; k<9; k++ ){
		if( min > v[k] ) min = v[k];
	}
	
	return min;
}

// cerca il tempo di collisione di un disco con le 9 repliche di un secondo
double amongst_copies_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s ){
	int p, q, k=0;
	double vec[ 100 ];
	
	for( p=-1; p<=1; p++){
		for( q=-1; q<=1; q++ ){
			vec[ k ] = coll_time_search( r1x, r1y, v1x, v1y, r2x + p, r2y + q, v2x, v2y, s );
			k++;
		}
	}
	
	/*
	k--;
	printf( "\n\n" );
	while( k>=0 ){
		printf( "%d___ %e ___\n", k, vec[k] );
		k--;
	}
	printf( "\n\n" );
	*/
	
	return v_min_search( vec );
}

// evoluzione del sistema
void evolve( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], double delta ){
	int k; // contatore del ciclo

	for( k = 0; k < N; k++ ){
		rx[ k ] = rx[ k ] + vx[ k ] * delta;
		ry[ k ] = ry[ k ] + vy[ k ] * delta;

		// proietto le posizioni dentro alla scatola
		rx[ k ] -= floor( rx[ k ] );
		ry[ k ] -= floor( ry[ k ] ); 
	}
}

// aggiornamento delle velocità
void speed_refresh( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], int p1, int p2 ){
	double vers_rx, vers_ry, dot_prod, sq_norm;

	// qualche calcolo preliminare	
	sq_norm = ( rx[ p2 ] - rx [ p1 ] ) * ( rx[ p2 ] -rx[ p1 ] ) + ( ry[ p2 ] -ry[ p1 ] ) * ( ry[ p2 ] -ry[ p1 ] );
	vers_rx = ( rx[ p2 ] - rx[ p1 ] ) / sqrt( sq_norm );
	vers_ry = ( ry[ p2 ] -ry[ p1 ] ) / sqrt( sq_norm );
	dot_prod = ( vx[ p1 ] - vx[ p2 ] ) * vers_rx + ( vy[ p1 ] -vy[ p2 ] ) * vers_ry;
	
	// assegno le nuvoe velocità
	vx[ p1 ] -= dot_prod * vers_rx;
	vy[ p1 ] -= dot_prod * vers_ry;
	vx[ p2 ] += dot_prod * vers_rx;
	vy[ p2 ] += dot_prod * vers_ry;
}
