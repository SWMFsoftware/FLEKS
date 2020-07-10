#ifndef _CONVERTER_H_
#define _CONVERTER_H_

#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_Box.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

class Converter{
private:
std::string dirIn; 
amrex::MultiFab mf; 
amrex::Box minBox; 
amrex::RealBox domainRange;
amrex::Geometry geom; 
std::string plot_string; 
amrex::Real rPlanet; 

int iCycle; 
amrex::Real time; 

amrex::Vector<std::string> varNames;

public:
Converter(const std::string& in){
    dirIn = in; 
} 
~Converter() = default; 

void read(); 
void convert();
void write(); 

};

#endif