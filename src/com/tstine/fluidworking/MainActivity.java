package com.tstine.fluidapp;

import android.app.Activity;
import android.os.Bundle;
import android.opengl.GLSurfaceView;
import android.os.Bundle;
import android.content.Context;

import android.view.MotionEvent;

import java.lang.Float;
import java.lang.Double;

import android.util.Log;


public class MainActivity extends Activity
{
	private GLSurfaceView gl_view;
	@Override
	public void onCreate(Bundle savedInstanceState)
	{
		super.onCreate(savedInstanceState);
		gl_view = new FluidSurface(this);
		setContentView(gl_view);
	}
}

class FluidSurface extends GLSurfaceView{
	FluidRenderer renderer;
	public FluidSurface(Context context){
		super(context);
		renderer = new FluidRenderer(context);
		setEGLContextClientVersion(2);
		setRenderer( renderer );
		//setRenderMode(GLSurfaceView.RENDERMODE_WHEN_DIRTY);
	}

	@Override
	public boolean onTouchEvent(MotionEvent e){
		float x = e.getX();
		float y = e.getY();
		int input = e.getAction();
		if( input == MotionEvent.ACTION_DOWN ){
			renderer.omy = y;
			renderer.omx = x;
			renderer.my = y;
			renderer.mx = x;

			renderer.touched = true;
		}
		if( input == MotionEvent.ACTION_UP ){
			renderer.touched = false;
		}
		if( input == MotionEvent.ACTION_MOVE ){
			renderer.my = y;
			renderer.mx = x;
		}
		/*		switch(e.getAction()){
		case MotionEvent.ACTION_MOVE:
			//Log.d("FluidRenderer", "size of float: " + Float.SIZE);
			//Log.d("FluidRenderer", "size of double: " + Double.SIZE);
			//Log.d( "FluidRenderer", "Moving! ( " + x + ", " + y + " )" );
			break;
		case MotionEvent.ACTION_DOWN:
			Log.d( "FluidRenderer", "Down! ( " + x + ", " + y + " )" );
			break;
		case MotionEvent.ACTION_UP:

			Log.d( "FluidRenderer", "Up! ( " + x + ", " + y + " )" );
			break;
			}*/


		return true;
	}
}

