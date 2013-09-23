package com.tstine.fluidapp;
  
import android.os.Bundle;

import android.app.Activity;

import javax.microedition.khronos.egl.EGLConfig;
import javax.microedition.khronos.opengles.GL10;

import android.opengl.GLSurfaceView;
import android.opengl.GLES20;
import android.opengl.GLU;
import android.opengl.Matrix;

import java.nio.IntBuffer;
import java.nio.FloatBuffer;
import java.nio.ByteBuffer;
import java.nio.ShortBuffer;
import java.nio.*;

import android.util.Log;

import java.io.*;
import java.lang.*;
import java.util.Calendar;
import java.util.Random;
import java.lang.Math;

import android.content.Context;
import android.content.res.AssetManager;

import java.lang.Exception;
import java.lang.StackTraceElement;


public class FluidRenderer implements GLSurfaceView.Renderer{
	public static final String TAG = "FluidRenderer";


	/** Fluid simulator variables
	 */

	public enum bnd_type{D_BND, U_BND, V_BND}
	public static final int N = 56;
	public static final int size = (N+2) * (N+2);
	public static float dt;
	public static float diffusivity;
	public static float viscosity;
	public static float force;
	public static float source_add;

	public static float death_rate;
 	
	public static float[] u;
	public static float[] v;
	public static float[] d;
	public static float[] u0;
	public static float[] v0;
	public static float[] d0;
	//new stuff
	public static float[] p;
	public static float[] p0;
	public static float[] RHS;
	public static float[] r;
	public static float[] fn;
	public static float[] gn;
	public static float[] Gx;
	public static float[] Gy;
	
	public static int Re;
	public static int max_iters;
	public static float tau;
	public static float w;
	public static float TOL;
	public static float dx;
	public static float dy;
	public static float gamma;
	/**End**/


	public Context ctx;
	public float mx=0.0f, my=0.0f, omx=0.0f, omy=0.0f, dmx = 0.0f, dmy = 0.0f;

	public static boolean touched = false;
	public static boolean add_density_source = false;
	public static boolean add_velocity = true;
	
	private Calendar calendar;
	private static int frames = 0;
	private static long T0 = 0;

	private float[] mvp_matrix = new float[16];
	private float[] projection_matrix = new float[16];
	private float[] modelview_matrix = new float[16];
	private float[] normal_matrix = new float[9];

	private int mvp_matrix_handle;

	public int m_width = 0;
	public int m_height = 0;

	private static float[] vertices;
  
	private static float[] colors;
	private static final float[] color = { 1.0f, 0.0f, 0.0f, 1.0f };
	private static FloatBuffer vert_buffer;
	private static FloatBuffer color_buffer;

	private static final int float_size = 4;
	private static final int short_size = 2;
	private static final int int_size = 4;

	private static final int num_indices = 3 * 2 * (N-1) * (N-1);
	private static int [] indices;
	private static IntBuffer indices_buffer;
	private static final int indices_attrib_len = 3;
	private static final int indices_size = num_indices *
		indices_attrib_len * int_size;

	private final static int num_vertices = N * N;
	private static float[] vertex_data;
	private static FloatBuffer vertex_buffer;

	private static int[] vbo_ids = new int[2];
	
	private static final int vertex_pos_len = 3;
	private static final int vertex_color_len = 3;
	//private static final int vertex_normal_len = 3;
	//private static final int vertex_texcoord_len = 2;

	private static final int vertex_pos_idx = 0;
	private static final int vertex_color_idx = 1;
	//private static final int vertex_normal_idx = 1;
	//private static final int vertex_texcoord_idx = 3;

	private static final int vertex_attrib_len = vertex_pos_len + vertex_color_len;
	//+vertex_normal_len + vertex_color_len + vertex_texcoord_len;
	
	private static final int vertex_pos_offset = 0;
	private static final int vertex_color_offset = 3;
	//private static final int vertex_normal_offset = 3;
	//private static final int vertex_texcoord_offset = 9;

	private static final int vertex_stride = vertex_attrib_len * float_size;
	private static final int vertex_array_len = vertex_attrib_len * num_vertices;
	private static final int vertex_array_size = vertex_attrib_len * num_vertices * float_size;
		
	private static int position_handle, color_handle, normal_handle, texcoord_handle;



	
	
	public FluidRenderer(Context context){
		this.ctx = context;
		Shader.ctx = context;
		try{
			vertex_data = new float[ vertex_array_len ];
			indices = new int[ num_indices ];
		}catch(Exception e){
			Log.e( TAG, (e.getStackTrace()).toString() );
		}catch( OutOfMemoryError e ){
			Log.e(TAG, "Out of memory");
		}
  
		vertex_buffer = ByteBuffer.allocateDirect(  vertex_array_size )
				.order( ByteOrder.nativeOrder() ).asFloatBuffer();
 
		indices_buffer = ByteBuffer.allocateDirect( indices_size )
		.order( ByteOrder.nativeOrder( ) ).asIntBuffer();
		
		set_default_values();
		allocate_memory();
	}

	private static void allocate_memory(){
		try{
			u = new float[size];
			v = new float[size];
			d = new float[size];
			u0 = new float[size];
			v0 = new float[size];
			d0 = new float[size];
			
			//new stuff
			p = new float[size];
			p0 = new float[size];
			Gx = new float[size];
			Gy = new float[size];
			fn = new float[size];
			gn = new float[size];
			RHS = new float[size];
			r = new float[size];
			
		}catch(Exception e){
			Log.e( TAG, (e.getStackTrace()).toString() );
		}catch( OutOfMemoryError e ){
			Log.e(TAG, "Out of memory");
		}
	}

	public static void set_default_values(){
		dt = .1f;
		diffusivity = .0000f;
		viscosity = 0.000f;
		force = 5.0f;
		source_add = 500.0f;
		death_rate = .1f;
		
		//new stuff
		max_iters = 20;
		tau = 0.5f;
		w = 1.7f;
		Re = 1000;
		TOL = .01f;
		dy = dx = 1.0f/(N+1);
		gamma = 0.9f;
	}

	@Override
	public void onSurfaceCreated(GL10 gl, EGLConfig config){

		GLES20.glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		GLES20.glGenBuffers( 2, vbo_ids, 0 ); 

		Shader.load();
		create_vertices();

		indices_buffer.put( indices ).position( 0 );
		vertex_buffer.put( vertex_data ).position( 0 );

		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, vbo_ids[1] );
		GLES20.glBufferData( GLES20.GL_ELEMENT_ARRAY_BUFFER, indices_size,
												 indices_buffer, GLES20.GL_STATIC_DRAW );

		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, vbo_ids[0] );
		GLES20.glBufferData( GLES20.GL_ARRAY_BUFFER, vertex_array_size,
												 vertex_buffer, GLES20.GL_DYNAMIC_DRAW );

		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, 0 );
		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, 0 );

		position_handle = GLES20.glGetAttribLocation( Shader.program_handle,
																									"vPosition");
		color_handle = GLES20.glGetAttribLocation( Shader.program_handle,
																							 "vColor");
		normal_handle = GLES20.glGetAttribLocation( Shader.program_handle,
																								"vNormal" );
		texcoord_handle = GLES20.glGetAttribLocation( Shader.program_handle,
																									"vTexcoord" );
		Log.d(TAG, "N= " + N );
	}
	
	@Override 
	public void onDrawFrame(GL10 gl){
		
		render(); 
		get_from_ui();
		step();
		calculate_fps(); 

	}

	public static int AT(int i, int j){return i+(N+2)*j;}

	private static void swap_d(){
		float[]temp = d;
		d = d0;
		d0 = temp;
	}
	private static void swap_u(){
		float[]temp = u;
		u = u0;
		u0 = temp;
	}
	
	private static void swap_v(){
		float[]temp = v;
		v = v0;
		v0 = temp;
	}

	public static void step(){
		step_density();
		step_velocity();
		kill_source();
	}

	private static void step_density(){
		add_source( d, d0 );

		swap_d( );
		diffuse( bnd_type.D_BND, d, d0, diffusivity );

		swap_d( );
		advect( bnd_type.D_BND, d, d0, u, v );
	}


	public static float d2x( float[] x, int i, int j ){
		return ( x[ AT(i+1,j)] - 2 * x[ AT(i,j) ] +
						 x[ AT(i-1,j) ] ) / (float) Math.pow(dx,2);
	}
	public static float d2y( float[] x, int i, int j ){
		return ( x[ AT(i,j+1) ] - 2*x[ AT(i,j) ]+
						 x[ AT(i,j-1) ] ) / (float)Math.pow(dy,2);
	}
	public static float du2dx( float[] x, int i, int j){
		return (1/dx * ((float)Math.pow( ( x[AT( i , j )] + x[AT( i+1 , j )])/2,2) -
										(float)Math.pow( ( x[AT( i-1 , j )] + x[AT( i , j )] )/2,2)) + 
						gamma/dx * (Math.abs( x[AT( i , j )] + x[AT( i+1 , j )])/2 *
												(x[AT( i , j )] - x[AT( i+1 , j )])/2 -
												Math.abs( x[AT( i-1 , j )] + x[AT( i , j )] )/2 *
												(x[AT( i-1 , j )] - x[AT( i , j )])/2));
	}

	public static float duvdy(float[] x, float[] y,
																		 int i, int j){
		return (1/dy *((y[AT( i , j )] + y[AT( i+1 , j )])/2 *
									 (x[AT( i , j )] + x[AT( i , j+1 )])/2 -
									 (y[AT( i , j-1 )] + y[AT( i+1 , j-1 )])/2 *
									 (x[AT( i , j-1 )] + x[AT( i , j )])/2) +
						gamma/dy *(Math.abs( y[AT( i , j )] + y[AT( i+1 , j )] )/2 *
											 (x[AT( i , j )] - x[AT( i , j+1 )])/2 -
											 Math.abs( y[AT( i , j-1 )] + y[AT( i , j )])/2 *
											 (x[AT( i , j-1 )] - x[AT( i , j )])/2));
	}

	public static float duvdx( float[] x, float[] y,
														 int i, int j){
		return (1/dx *((y[AT( i , j )] + y[AT( i+1 , j )])/2 *
									 (x[AT( i , j )] + x[AT( i , j+1 )])/2 -
									 (y[AT( i-1 , j )] + y[AT( i , j )])/2 *
									 (x[AT( i-1 , j )] + x[AT( i-1 , j+1 )])/2) +
						gamma/dx *(Math.abs( x[AT(i,j)] + x[AT(i,j+1)] )/2 *
											 (y[AT(i,j)] - y[AT(i+1,j)])/2 -
											 Math.abs(x[AT(i=1,j)] + x[AT(i-1,j+1)])/2 *
											 (y[AT(i-1,j)] - y[AT(i,j)])/2));
	}

	public static float dv2dy( float[] y, int i, int j){
		return ( 1/dy * ( (float)Math.pow( ( y[ AT(i,j) ] + y[ AT(i,j+1) ] )/2, 2 ) -
											(float)Math.pow( ( y[ AT(i,j-1) ] + y[ AT(i,j) ] )/2, 2) ) + 
						 gamma/dy * ( Math.abs( y[ AT(i,j) ] + y[ AT(i,j+1) ] )/2 *
													( y[ AT(i,j )] - y[ AT(i,j+1) ] )/2 -
													Math.abs( y[ AT(i,j-1) ] + y[ AT(i,j) ] )/2 *
													( y[ AT(i,j-1) ] - y[ AT(i,j) ] )/2 ) );
	}

	public static int eE(int i){return i<N? 1:0;}
	public static int eW(int i){return i==0? 0:1;}
	public static int eN(int j){return j<N? 1:0;}
	public static int eS(int j){return j==0? 0:1;}


	private static void step_velocity(){

		/*
		int i,j;
		int imax = N, jmax = N;

		for( j = 1; j <= jmax; j++ ){
			i=j;
			v[ AT( 0,j ) ] = -v[ AT( 1,j ) ];
			v[ AT( imax+1,j ) ] = -v[ AT( imax,j ) ];
			u[ AT( i,0 ) ] = -u[ AT( i,1 ) ];
			u[ AT( i,jmax+1 ) ] = -u[ AT( i,jmax ) ];
		}

		for( i=1; i <= imax; i++ ){
			for( j=1; j <= jmax; j++ ){
				fn[ AT(i,j) ] = u[ AT(i,j) ] + dt * ( 1/Re * ( d2x(u,i,j) +
																											 d2y(u,i,j) ) -
																							du2dx(u,i,j) - duvdy(u,v,i,j) +
																							Gx[ AT(i,j) ] );

				gn[ AT(i,j) ] = v[ AT(i,j) ] + dt * ( 1/Re * ( d2x(v,i,j) +
																											 d2y(v,i,j) ) -
																							duvdx(u,v,i,j) - dv2dy(v,i,j) +
																							Gy[ AT(i,j) ] );
			}
		}

		for( i=1; i <=imax; i++ ){
			j = i;
			fn[ AT( 0,j ) ] = u[ AT( 0,j ) ];
			fn[ AT( imax,j ) ] = u[ AT( imax, j ) ];
			gn[ AT( i,0 ) ] = v[ AT( i,0 ) ];
			gn[ AT( i, jmax ) ] = v[ AT( i, jmax ) ];
		}

		for( i=1; i<=imax; i++ ){
			for( j=1; j<=jmax; j++ ){
				RHS[ AT( i,j ) ] = 1/dt * ( (fn[ AT( i,j ) ] -
																			fn[ AT ( i-1,j ) ])/dx +
																		(gn[ AT( i,j ) ] -
																			gn[ AT( i,j-1 ) ])/dy );

			}
		}

		int it = 0;
		while( it < max_iters ){
			for( i=1; i<=imax; i++ ){
				for( j=1; j<=jmax; j++ ){
					p[ AT( i,j ) ] = (1-w) * p0[ AT( i,j ) ]+
						w/( ( eE(i) + eW(i) ) / (float)Math.pow( dx,2 )
								+ ( eN(j)+eS(j) ) / (float)Math.pow( dy,2 ) ) *
						( (eE(i) * p0[ AT( i+1,j ) ] + eW(i)
							 * p[ AT( i-1,j ) ] ) / (float)Math.pow( dx,2 ) +
							( eN(j) * p0[ AT( i,j+1 ) ] + eS(j)
								* p[ AT( i,j-1 ) ] ) / (float)Math.pow( dy,2 )
							- RHS[ AT( i,j ) ] );

							r[ AT( i,j ) ] = ( eE(i) *
														 ( p0[ AT( i+1,j ) ] - p0[ AT( i,j ) ] ) -
														 eW(i) * ( p0[ AT( i,j ) ] - p0[ AT( i-1,j ) ] ) )/Math.pow(dx,2) +
						( eN(j) * ( p0[ AT( i,j+1 ) ] - p0[ AT( i,j ) ] ) -
							eS(j) * ( p0[ AT( i,j ) ] - p0[ AT( i,j-1 ) ] ) ) / Math.pow(dy,2)
							- RHS[ AT( i,j ) ];
				}
			}
			it++;
		}

		float u_max = -1.0e10f, v_max= -1.0e10f;
		for( i=1; i<=imax; i++ ){
			for( j=1; j<=jmax; j++ ){
				u[ AT(i,j) ] = fn[ AT(i,j) ] - dt/dx *( p[ AT(i+1,j) ] - p[ AT(i,j) ] );
				v[ AT(i,j) ] = gn[ AT(i,j) ] - dt/dy *( p[ AT(i,j+1) ] - p[ AT(i,j) ] );
				if( u[ AT(i,j) ] > u_max ) u_max = u[ AT(i,j) ];
				if( v[ AT(i,j) ] > v_max ) v_max = v[ AT(i,j) ];
			}
		}

		float a = ( Re/2 )/( 1/(float)Math.pow(dx,2) + 1/(float)Math.pow(dy,2) );
		float b = dx / Math.abs( u_max );
		float c = dy / Math.abs(v_max);
		dt = tau * Math.min( Math.min( a,b ),c );
		*/

		add_source( u,u0 );
		add_source( v,v0 );

		swap_u();
		diffuse( bnd_type.U_BND, u, u0, viscosity );

		swap_v( );
		diffuse( bnd_type.V_BND, v, v0, viscosity );

		project( );
		
		swap_u( );
		swap_v( );

		advect( bnd_type.U_BND, u, u0, u0, v0 );
		advect( bnd_type.V_BND, v, v0, u0, v0 );

		project( );

	}
	
	private static void add_source(float[] x, float[] src){
		int i = 0;
		for( i=0; i < size; i++ ){
			x[i] += dt * src[i];
		}
		
	}

	private static void kill_source(){
		int i = 0;
		for( i=0; i< size; i++ ){
			d[i] = Math.max( 0.0f, ( d[i] - dt*death_rate ) );
		}
	}
	
	private static void diffuse( bnd_type bnd,
															 float[]x, float[] x0, float diff ){
		float num = dt * diff * N * N;
		linear_solver(bnd, x, x0, num, 1+4*num);
		
	}

	private static void linear_solver(bnd_type bnd, float[]x, float[]x0,
																		float num, float denom){
		int i=0,j=0,k=0;
		float avg_chg = 0.0f, val = 0.0f, max_change = -1.0f, change = -1.0f;;

		for(k=0; k<4; k++){
			for(i=1; i<=N; i++){
				for(j=1; j<=N; j++){
					val = (x0[AT(i,j)] + num *(x[AT(i-1, j)] + x[AT(i+1,j)]+
																					x[AT(i, j-1)] + x[AT(i,j+1)] ))/denom;
					change = Math.abs( val - x[AT(i,j)] );
					if( change > max_change ){
						max_change = change;
					}
					x[AT(i, j)] = val;
				}
			}
			set_bnd(bnd, x);
			max_change = -1.0f;
		}


	}

	private static void advect(bnd_type bnd, 
														 float[] x, float[] x0,
														 float[] u, float[] v){
		int i, j, i1, j1, i0, j0;
		float dt0, x_pos, y_pos, s0, s1, t0, t1;
		dt0 = dt * N;
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){

				x_pos = i - dt0 * u[AT(i,j)];
				y_pos = j - dt0 * v[AT(i,j)];
				
				if( x_pos < 0.5f )
					x_pos = 0.5f;
				else if( x_pos > N+0.5f )
					x_pos = N + 0.5f;
				i0=(int)x_pos;
				i1 = i0+1;
				
				if( y_pos < 0.5f )
					y_pos = 0.5f;
				if( y_pos > N+0.5f )
					y_pos = N + 0.5f;
				j0=(int)y_pos;
				j1 = j0+1;
				
				s1 = x_pos - i0;
				s0 = 1 - s1;
				t1 = y_pos - j0;
				t0 = 1 - t1;
				
				x[AT(i,j)] = s0*( t0 * x0[AT( i0,j0 )] + t1 * x0[AT( i0,j1 )] ) + 
					s1 * ( t0 * x0[AT( i1,j0 )] + t1 * x0[AT( i1,j1 )] );
			}
		}

		set_bnd( bnd,x );
	}

	private static void project(){
		int i=0, j=0;
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){
				//p=u0, div=v0
				v0[ AT( i,j )] = -0.5f *( u[AT( i+1,j )] - u[AT( i-1,j )] +
																	v[AT( i,j+1 )] - v[AT( i,j-1 )] )/N;
				u0[AT( i,j )] = 0;
			}
		}
		set_bnd( bnd_type.D_BND, u0 );
		set_bnd( bnd_type.D_BND, v0 );
		linear_solver( bnd_type.D_BND, u0, v0, 1, 4 );
		
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){
				u[AT( i,j )] -= 0.5f * N * ( u0[AT( i+1,j )]-u0[AT( i-1,j )] );
				v[AT( i,j )] -= 0.5f * N * ( u0[AT( i,j+1 )]-u0[AT( i,j-1 )] );
			}
		}
		set_bnd(bnd_type.U_BND, u);
		set_bnd(bnd_type.V_BND, v);
	}

	private static void set_bnd(bnd_type type, float[] x){
		int i=0;
		for(i=1; i<=N; i++){
			switch( type ){
			case U_BND:
				x[AT( 0,i )]   = -x[AT( 1,i )];
				x[AT( N+1,i )] = -x[AT( N,i )];
				x[AT( i,1 )]   = 0.0f;
				x[AT( i,N )]   = 0.0f;
				break;
			case V_BND:
				x[AT( 0,i )]   = 0.0f;
				x[AT( N+1,i )] = 0.0f;
				x[AT( i,0 )]   = -x[AT( i,1 )];
				x[AT( i,N+1 )] = -x[AT( i,N )];
				break;
			case D_BND:
				x[AT( 0,i )]   = x[AT( 1,i )];
				x[AT( N+1,i )] = x[AT( N,i )];
				x[AT( i,0 )]   = x[AT( i,1 )];
				x[AT( i,N+1 )] = x[AT( i,N )];
				break;
			}
		}
		x[AT( 0,0 )]     = 0.5f * ( x[AT( 1,0 )] + x[AT( 0,1 )] );
		x[AT( 0,N+1 )]   = 0.5f * ( x[AT( 1,N+1 )] + x[AT( 0,N )] );
		x[AT( N+1,0 )]   = 0.5f * ( x[AT( N,0 )] + x[AT( N+1,1 )] );
		x[AT( N+1,N+1 )] = 0.5f * ( x[AT( N,N+1 )] + x[AT( N+1,N )] );
	}


 
	@Override
	public void onSurfaceChanged(GL10 gl, int width, int height){
		m_width = width;
		m_height = height;
		GLES20.glViewport( 0, 0, width, height );
		float aspect = (float)width/height;
		Matrix.setIdentityM(projection_matrix, 0);
		Matrix.orthoM( projection_matrix,
									 0,
									 -1, 1,
									 -1, 1,
									 0.1f, 100.0f );
	}


	public void get_from_ui(){

		int i;
		for(i=0; i< size; i++){
			u0[i] = v0[i] = d0[i] = 0.0f;
		}

		float delta = 5.0f;
		float dist = ( float )Math.sqrt( Math.pow( mx - omx,2 ) + Math.pow( my - omy, 2) );

		if(touched && dist > delta){
			int x_cell = (int)( ( omx/(float)m_width )*N+1 );
			int y_cell = (int) ( ( (m_height-omy )/(float)m_height )*N+1 );

			x_cell = clamp( x_cell, 1, N );
			y_cell = clamp( y_cell, 1, N );
			if(x_cell >= 1 || x_cell <= N || y_cell >= 1 || y_cell <= N) {
				d0[AT( x_cell,y_cell )] = source_add;
				u0[AT( x_cell,y_cell )] = force * ( mx - omx );
				v0[AT( x_cell,y_cell )] = force * -( my - omy );
			}
		}
		omx = mx;
		omy = my;
	}

	public void render(){
		GLES20.glClear(GLES20.GL_COLOR_BUFFER_BIT);
		GLES20.glUseProgram( Shader.program_handle );

		set_mvp_matrix();

		//render_triangle();
		render_simulator();
	}

	public void set_mvp_matrix(){

		Matrix.setIdentityM(modelview_matrix, 0);

		Matrix.setLookAtM(modelview_matrix, 0,
											0, 0, 3,
											0, 0, 0,
											0, 1, 0);

		Matrix.multiplyMM(mvp_matrix, 0, projection_matrix, 0, modelview_matrix, 0);

		mvp_matrix_handle = GLES20.glGetUniformLocation( Shader.program_handle, "mvp_matrix" );
		FluidRenderer.checkGlError( "glGetUniformLocation" );

		GLES20.glUniformMatrix4fv( mvp_matrix_handle, 1, false, mvp_matrix, 0 );
		FluidRenderer.checkGlError( "glUniformMatrix4fv" );		
	}

	public void render_simulator(){
		set_colors();
		vertex_buffer.put( vertex_data ).position( 0 );

		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, vbo_ids[1] );
		//GLES20.glBufferData( GLES20.GL_ELEMENT_ARRAY_BUFFER, indices_size,
												 //indices_buffer, GLES20.GL_STATIC_DRAW );

		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, vbo_ids[0] );
		//GLES20.glBufferSubData( GLES20.GL_ARRAY_BUFFER, 0,
		//vertex_array_size, vertex_buffer );
		GLES20.glBufferSubData( GLES20.GL_ARRAY_BUFFER, 0, 
														vertex_array_size, vertex_buffer );


		GLES20.glEnableVertexAttribArray( position_handle );
		GLES20.glEnableVertexAttribArray( color_handle );
		//GLES20.glEnableVertexAttribArray( normal_handle );
		//GLES20.glEnableVertexAttribArray( texcoord_handle );

		//position information
		GLES20.glVertexAttribPointer( position_handle, vertex_pos_len,
																	GLES20.GL_FLOAT, false,
																	vertex_stride,
																	vertex_pos_offset * float_size );
		//color information
		GLES20.glVertexAttribPointer( color_handle, vertex_color_len,
																	GLES20.GL_FLOAT, false,
																	vertex_stride,
																	vertex_color_offset * float_size );


		/*		//normal information
		GLES20.glVertexAttribPointer( normal_handle, vertex_normal_len,
																	GLES20.GL_FLOAT, false,
																	vertex_stride,
																	vertex_normal_offset * float_size );

		//texcoord information
		GLES20.glVertexAttribPointer( texcoord_handle, vertex_texcoord_len,
																	GLES20.GL_FLOAT, false,
																	vertex_stride,
																	vertex_texcoord_offset * float_size );
		*/

		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, vbo_ids[1] );
		GLES20.glDrawElements( GLES20.GL_TRIANGLES,
													 num_indices,
													 GLES20.GL_UNSIGNED_INT,
													 0 );

		GLES20.glDisableVertexAttribArray( position_handle );
		GLES20.glDisableVertexAttribArray( color_handle );
		//GLES20.glDisableVertexAttribArray( normal_handle );
		//GLES20.glDisableVertexAttribArray( texcoord_handle );

 
		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, 0 );
		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, 0 );


	}
	public void set_colors(){
		int i,j,c_idx=0, idx=0;
		float dens;
		Random rand = new Random();
		for(j=1; j< N+1; j++){
			for(i=1; i< N+1; i++){
				dens = d[AT(i,j)];
				c_idx = ( idx * vertex_attrib_len ) + vertex_color_offset;
				
				vertex_data[ c_idx ] = dens;
				vertex_data[ c_idx+1 ] = dens;
				vertex_data[ c_idx+2 ] = dens;
				idx++;
			}
		}
	}

	public void render_triangle(){
		float[] tri_vertices = {-1.0f, -1.0f, 0.0f,
														1.0f,  0.0f, 0.0f,
														-1.0f, 1.0f, 0.0f};
		float[] tri_colors = {1.0f, 0.0f, 0.0f, 1.0f,
													0.0f, 1.0f, 0.0f, 1.0f,
													0.0f, 0.0f, 1.0f, 1.0f};
		int coords_per_vertex = 3;
		int colors_per_vertex = 4;
		ByteBuffer v_bb = ByteBuffer.allocateDirect( tri_vertices.length * 4 );
		ByteBuffer c_bb = ByteBuffer.allocateDirect( tri_colors.length * 4 );
		v_bb.order(ByteOrder.nativeOrder());
		c_bb.order(ByteOrder.nativeOrder());
		FloatBuffer vert_buffer = v_bb.asFloatBuffer();
		FloatBuffer color_buffer = c_bb.asFloatBuffer();
		vert_buffer.put(tri_vertices);
		color_buffer.put(tri_colors);
		vert_buffer.position(0);
		color_buffer.position(0);

		int position_handle = GLES20.glGetAttribLocation(Shader.program_handle, "vPosition");
		int color_handle = GLES20.glGetAttribLocation(Shader.program_handle, "vColor");

		GLES20.glEnableVertexAttribArray( position_handle );
		GLES20.glEnableVertexAttribArray( color_handle );
		GLES20.glVertexAttribPointer( position_handle, coords_per_vertex,
																 GLES20.GL_FLOAT, false,
																	coords_per_vertex *  4, vert_buffer);

		GLES20.glVertexAttribPointer( color_handle, colors_per_vertex,
																GLES20.GL_FLOAT, false,
																	colors_per_vertex * 4, color_buffer );

		GLES20.glDrawArrays( GLES20.GL_TRIANGLES, 0, 3 );
		GLES20.glDisableVertexAttribArray( position_handle );
		GLES20.glDisableVertexAttribArray( color_handle );
	}

	public void calculate_fps(){
		frames++;
		calendar = Calendar.getInstance();
		long t = calendar.getTimeInMillis();
		if(t - T0 >= 5000){
			float seconds = (float)(t-T0)/1000.0f;
			float fps = frames/seconds;
			Log.d(TAG, frames + " frame in "+ seconds + " seconds = " +
						fps + " FPS\n");
			T0 = t;
			frames = 0;
		}
	}

	private static int v_idx_at(int i, int j){ return  (i+N*j); }

	private static void create_vertices(){
		int i,j;
		int vert_idx=0, norm_idx=0, tex_idx = 0, i_idx=0, idx=0;
		float h = 2.0f/(N-1);
		float x,y,d;
		for(j=0; j< N; j++){
			y = -1 + (j * h);
			for(i=0; i< N; i++){
				x = -1 + (i * h);
				//Log.d(TAG, "x,y= ( " + x + ", " + y + " )");

				vert_idx = ( idx * vertex_attrib_len ) + vertex_pos_offset; 
				//norm_idx = ( idx * vertex_attrib_len ) + vertex_normal_offset; 
				//tex_idx = ( idx * vertex_attrib_len ) + vertex_texcoord_offset; 

				vertex_data[ vert_idx   ] = x;
				vertex_data[ vert_idx+1 ] = y;
				vertex_data[ vert_idx+2 ] = 0.0f;

				//vertex_data[ norm_idx ] = 0.0f;
				//vertex_data[ norm_idx+1 ] = 0.0f;
				//vertex_data[ norm_idx+2 ] = 1.0f;

				//vertex_data[ tex_idx ] = x;
				//vertex_data[ tex_idx+1 ] = y;
				
				idx++;
			}
		}

		for(j=0; j < N-1; j++){
			for(i=0; i < N-1; i++){
				indices[i_idx++] = v_idx_at( i,j );
				indices[i_idx++] = v_idx_at( i+1,j );
				indices[i_idx++] = v_idx_at( i+1,j+1 );
 
				indices[i_idx++] = v_idx_at( i,j );
				indices[i_idx++] = v_idx_at( i+1,j+1 );
				indices[i_idx++] = v_idx_at( i,j+1 );
			}
		}
	}

    public static void checkGlError(String glOperation) {
        int error;
        while ((error = GLES20.glGetError()) != GLES20.GL_NO_ERROR) {
            Log.e(TAG, glOperation + ": glError " + error);
            throw new RuntimeException(glOperation + ": glError " + error);
        }
    }

	public static int clamp(int value, int min, int max){
		value = value > max ? max : value;
		value = value < min ? 1 : value;
		return value;
	}


}

/*
	 Vertex buffer object class, contains all of the code
	 for drawing a vertex buffer.  This is mostly specalized
	 to this particular usage of a VBO and would need
	 to be modified to be more generic
*/

/*
class VertexBufferObject{

	private static class Attrib{
		public static int POSITION= 0;
		public static int NORMAL= 1;
		public static int COLOR= 2;
		public static int TEXCOORD= 3;
		public static final int SIZE = 4;
	}

	private static final String TAG = "FluidRenderer";

	private final int SIZEOF_FLOAT = 4;
	private final int SIZEOF_INT = 4;

	private int[] mLength;
	private int[] mOffset;

	private int mVertex_length;
	private int mVertex_stride;

	private float [] mVertex_data;
	private int [] mIndices_data;

	private FloatBuffer mVertex_buffer;
	private IntBuffer mIndices_buffer;

	private int mN;

	private int mNum_vertices;
	private int mNum_indices;

	private int[] mVbo_ids;

	private int mPosition_handle;
	
	public VertexBufferObject( int dimension ){
		this.mN = dimension;

		this.mNum_vertices = mN * mN;
		this.mNum_indices = 3 * 2 * (mN-1) * (mN -1);

		this.mLength = new int[ Attrib.SIZE ];
		this.mOffset = new int[ Attrib.SIZE ];
		this.mVbo_ids = new int[2];

		this.mLength[ Attrib.POSITION ] = 3;
		this.mLength[ Attrib.NORMAL ] = 3;
		this.mLength[ Attrib.COLOR ] = 3;
		this.mLength[ Attrib.TEXCOORD ] = 2;

		this.mVertex_length = this.mLength[ Attrib.POSITION ];

		this.mOffset[ Attrib.POSITION ] = 0;
		this.mOffset[ Attrib.NORMAL ] = 3;
		this.mOffset[ Attrib.COLOR ] = 6;
		this.mOffset[ Attrib.TEXCOORD ] = 9;
		
		this.mVertex_stride = mVertex_length * SIZEOF_FLOAT;

		try{
			mVertex_data = new float[ mVertex_length * mNum_vertices ];
			mIndices_data = new int[ mNum_indices ];
			mVbo_ids = new int[2];
		}catch(Exception e){
			Log.e( TAG, (e.getStackTrace()).toString() );
		}catch( OutOfMemoryError e ){
			Log.e(TAG, "Out of memory");
		}
  
		mVertex_buffer = ByteBuffer.allocateDirect(  mVertex_length * mNum_vertices * SIZEOF_FLOAT )
				.order( ByteOrder.nativeOrder() ).asFloatBuffer();
 
		mIndices_buffer = ByteBuffer.allocateDirect( mNum_indices * SIZEOF_INT )
		.order( ByteOrder.nativeOrder( ) ).asIntBuffer();

		GLES20.glGenBuffers( 2, mVbo_ids, 0 );

		this.create_vertices();
		this.create_indices();

		mIndices_buffer.put( mIndices_data ).position( 0 );
		mVertex_buffer.put( mVertex_data ).position( 0 );
		
		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, mVbo_ids[1] );
		GLES20.glBufferData( GLES20.GL_ELEMENT_ARRAY_BUFFER, mNum_indices * SIZEOF_INT,
												 mIndices_buffer, GLES20.GL_STATIC_DRAW );

		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, mVbo_ids[0] );
		GLES20.glBufferData( GLES20.GL_ARRAY_BUFFER, mVertex_length * SIZEOF_FLOAT,
												 mVertex_buffer, GLES20.GL_STATIC_DRAW );

		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, 0 );
		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, 0 );

		if(!Shader.loaded){
			throw new IllegalStateException("Shader not loaded!"+
																			"Call Shader.load() before this constructor");
		}

		mPosition_handle = GLES20.glGetAttribLocation( Shader.program_handle,
																									"vPosition");
	}

	private void create_vertices(){
		int i,j;
		int vert_idx=0, idx=0;
		float h = 2.0f/(mN-1);
		float x,y,d;
		for(j=0; j< mN; j++){
			y = -1 + (j * h);
			for(i=0; i< mN; i++){
				x = -1 + (i * h);

				vert_idx = ( idx * mVertex_length ) + mOffset[ Attrib.POSITION ];
				
				mVertex_data[ vert_idx   ] = x;
				mVertex_data[ vert_idx+1 ] = y;
				mVertex_data[ vert_idx+2 ] = 0.0f;
				
				idx++;
			}
		}
	}

	private int indices_at(int i, int j){return  (i+mN*j); }
	private void create_indices(){
		int i_idx=0, i, j;
		for(j=0; j < mN-1; j++){
			for(i=0; i < mN-1; i++){
				mIndices_data[i_idx++] = indices_at( i,j );
				mIndices_data[i_idx++] = indices_at( i+1,j );
				mIndices_data[i_idx++] = indices_at( i+1,j+1 );
 
				mIndices_data[i_idx++] = indices_at( i,j );
				mIndices_data[i_idx++] = indices_at( i+1,j+1 );
				mIndices_data[i_idx++] = indices_at( i,j+1 );
			}
		}
	}
	
	public void draw(){
		GLES20.glBindBuffer( GLES20.GL_ELEMENT_ARRAY_BUFFER, mVbo_ids[1] );
		GLES20.glBindBuffer( GLES20.GL_ARRAY_BUFFER, mVbo_ids[0] );

		GLES20.glEnableVertexAttribArray( mPosition_handle );

		//position information
		GLES20.glVertexAttribPointer( mPosition_handle, 3,
																	GLES20.GL_FLOAT, false,
*/

/**
	 Shader class, everything is prepared and loaded
	 here for the shader
 */
class Shader{
	public static int program_handle;
	public static Context ctx;
	public static boolean loaded=false;

	public static void load(){
		Shader.load("vertex_shader.glsl", "fragment_shader.glsl");
		loaded = true;
	}
	public static void load(String vertex_filename, String fragment_filename){
		try{
			Shader.vertex_code = load_shader_code(vertex_filename);
		}catch (IOException e){
			Log.e(TAG, "ERROR OPENING SHADER CODE FILE");
		}
		try{
			Shader.fragment_code = load_shader_code(fragment_filename);
		}catch(IOException e){
			Log.e(TAG, "ERROR OPENING FRAGMENT CODE FILE");
		}
		Shader.vertex_handle = load_shader(Shader.vertex_code,
																			 GLES20.GL_VERTEX_SHADER);
		Shader.fragment_handle = load_shader(Shader.fragment_code,
																				 GLES20.GL_FRAGMENT_SHADER);
		Shader.load_program();
	}

	private static final String TAG = "FluidRenderer";
	private static int vertex_handle;
	private static String vertex_code;

	private static String fragment_code;
	private static int fragment_handle;

	public static String load_shader_code(String filename) throws IOException{
		BufferedReader reader =
			new BufferedReader(new InputStreamReader(ctx.getAssets().open(filename)));
		StringBuffer sb = new StringBuffer();
		String line = null;
		while((line = reader.readLine()) != null){
			sb.append(line).append("\n");
		}
		reader.close();
		return sb.toString();
	}
	
	private static void load_program() throws RuntimeException {
		Shader.program_handle = GLES20.glCreateProgram();
		GLES20.glAttachShader( Shader.program_handle, Shader.vertex_handle );
		GLES20.glAttachShader( Shader.program_handle, Shader.fragment_handle );
		GLES20.glLinkProgram( Shader.program_handle );
		int[] compiled = new int[1];
		GLES20.glGetProgramiv( Shader.program_handle, 
													 GLES20.GL_LINK_STATUS,
													 compiled,
													 0 );
		while( compiled[0] == GLES20.GL_FALSE ){
			Log.d(TAG, "Linking shader failed!");
			String program_failure_info = GLES20.glGetProgramInfoLog( Shader.program_handle);
			throw new RuntimeException( "linking failed " );
		}
	}

	private static int load_shader(String shader_code, int type) throws RuntimeException{
		int shader = GLES20.glCreateShader( type );
		GLES20.glShaderSource( shader, shader_code );
		GLES20.glCompileShader( shader );
		int[] compiled = { 10 };
		GLES20.glGetShaderiv( shader, GLES20.GL_SHADER_SOURCE_LENGTH, compiled, 0 );
		while( compiled[0] == GLES20.GL_FALSE ){
			Log.d(TAG, "Could not compile shader!");
			String shader_failure_info = GLES20.glGetShaderInfoLog( shader );
			Log.d(TAG, "failure info: " + shader_failure_info );
			GLES20.glDeleteShader( shader );
			throw new RuntimeException( "compillation failed " );
		}
		return shader;
	}
}