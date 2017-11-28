#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom2.h>
#include <THStack.h>
#include "ProcessDir.h"
#include <TVector3.h>
#include <TMath.h>
#include <time.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include <ctime>
#include <chrono>
#define PI 3.14159265
#define target_radius 0.318

using namespace std;
void PrimaryProton(Float_t& vec1_x, Float_t& vec1_y, Float_t& vec1_z, Float_t p1_x, Float_t p1_y, Float_t p1_z)
{
    TVector3 centre_xy(0.16,0.21,0.0);
    TVector3 p1_xy(p1_x,p1_y,0.0);
    TVector3 vec1(vec1_x,vec1_y,vec1_z);
    Double_t theta = acos(vec1(2)/vec1.Mag());
    Double_t det = tan(theta);
    TVector3 vec1_xy(vec1_x,vec1_y,0.0);
    vec1_xy = (1.0/vec1_xy.Mag()) * vec1_xy;
    if ((p1_xy+det*vec1_xy-centre_xy).Mag() <= target_radius)
    {
        TVector3 k(0.0,0.0,1.0);
        vec1 = det*vec1_xy+k;
    }
    else
    {
        //setting up quadratic equation coefficients
        Double_t a = 1.0;
        Double_t b = 2.0*(p1_xy-centre_xy).Dot(vec1_xy);
        Double_t c = (p1_xy-centre_xy).Mag2()-pow(target_radius,2.0);
        //selecting the positive distance (along the momentum x-y projection)
        Double_t d_xy = (-b + pow(pow(b,2.0) + 4*a*c,0.5)) / (2*a);
        TVector3 k(0.0,0.0,d_xy/det);
        vec1 = d_xy*vec1_xy+k;
    }
    vec1_x = vec1(0);
    vec1_y = vec1(1);
    vec1_z = vec1(2);
    return;
}

void ExitingHadron(Float_t& vec2_x, Float_t& vec2_y, Float_t& vec2_z, Float_t& p2_x, Float_t& p2_y, Float_t& p2_z)
{
    TVector3 centre_xy(0.16,0.21,0.0);
    TVector3 p2_xy(p2_x,p2_y,0.0);
    TVector3 vec2(vec2_x,vec2_y,vec2_z);
    Double_t theta = acos(vec2(2)/vec2.Mag());
    Double_t det = p2_z*tan(theta);
    TVector3 vec2_xy(vec2_x,vec2_y,0.0);
    vec2_xy = (1.0/vec2_xy.Mag()) * vec2_xy;
    if ((p2_xy-det*vec2_xy-centre_xy).Mag() <= target_radius)
    {
        TVector3 k(0.0,0.0,p2_z);
        vec2 = det*vec2_xy+k;   
    }
    else
    {  
        //setting up quadratic equation coefficients 
        Double_t a = 1.0;
        Double_t b = 2.0*(p2_xy-centre_xy).Dot(vec2_xy);
        Double_t c = (p2_xy-centre_xy).Mag2()-pow(target_radius,2.0);
        //selecting the negative distance (backwards wrt the momentum x-y projection)
        Double_t d_xy = (-b - pow(pow(b,2.0) + 4*a*c,0.5)) / (2*a);
        TVector3 k(0.0,0.0,p2_z-d_xy/det);
        vec2 = d_xy*vec2_xy+k;
    }
    p2_x += -1.0*vec2(0);
    p2_y += -1.0*vec2(1);
    p2_z += -1.0*vec2(2);
    vec2_x = vec2(0);
    vec2_y = vec2(1);
    vec2_z = vec2(2);
    return;
}

Double_t DistanceAtPointOfClosestApproach_Lines(Float_t vec1_x, Float_t vec1_y, Float_t vec1_z,
                                                       Float_t p1_x, Float_t p1_y, Float_t p1_z,
                                                       Float_t vec2_x, Float_t vec2_y, Float_t vec2_z,
                                                       Float_t p2_x, Float_t p2_y, Float_t p2_z,
                                                       Double_t& z)
{
    TVector3 v1(vec1_x,vec1_y,vec1_z);
    TVector3 v2(vec2_x,vec2_y,vec2_z);
    TVector3 point1(p1_x,p1_y,p1_z);
    TVector3 point2(p2_x,p2_y,p2_z);
    TVector3 dist;

    Double_t a = v1.Mag2();
    Double_t b = v1.Dot(v2);
    Double_t c = v2.Mag2();
    Double_t d = (point1 - point2).Dot(v1);
    Double_t e = (point1 - point2).Dot(v2);
    Double_t det = a*c - pow(b,2.0);
    Double_t sN, sD, tN, tD;
    sD = det;
    tD = det;
    if (fabs(det) < pow(10.0,-8.0))
    {   
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {   
        sN = (b * e - c * d);
        tN = (a * e - b * d);
    }
    if (fabs(sN) < pow(10.0,-8.0)) sN = 0.0;
    if (fabs(tN) < pow(10.0,-8.0)) tN = 0.0;
    dist = (point1 - point2) + (sN/sD)*v1 - (tN/tD)*v2;
    z = (point1 + (sN/sD)*v1)(2);
    return dist.Mag();   
}

Double_t DistanceAtPointOfClosestApproach_LineSegments(Float_t vec1_x, Float_t vec1_y, Float_t vec1_z,
                                                       Float_t p1_x, Float_t p1_y, Float_t p1_z, 
                                                       Float_t vec2_x, Float_t vec2_y, Float_t vec2_z,
                                                       Float_t p2_x, Float_t p2_y, Float_t p2_z,
                                                       Double_t& z)
{   
    //should call PrimaryProton and ExitingHadron functions to modify the vectors and points before calculating the distance of closest approach
    PrimaryProton(vec1_x,vec1_y,vec1_z,p1_x,p1_y,p1_z); 
    ExitingHadron(vec2_x,vec2_y,vec2_z,p2_x,p2_y,p2_z);
    
    TVector3 v1(vec1_x,vec1_y,vec1_z);
    TVector3 v2(vec2_x,vec2_y,vec2_z);
    TVector3 point1(p1_x,p1_y,p1_z);
    TVector3 point2(p2_x,p2_y,p2_z);
    TVector3 dist;
    //magnitude squared of vector parellel to incoming particle direction
    Double_t a = v1.Mag2();
    //dot product of vectors parallel to incoming and outgoing particle directions
    Double_t b = v1.Dot(v2);
    //magnitude squared of vector parellel to outgoing particle direction
    Double_t c = v2.Mag2();
    Double_t d = (point1 - point2).Dot(v1);
    Double_t e = (point1 - point2).Dot(v2);
    //check if the incoming and outgoing directions are non-parallel
    Double_t det = a*c - pow(b,2.0);
    Double_t sN, sD, tN, tD;
    sD = det;
    tD = det;
    if (fabs(det) < pow(10.0,-8.0))
    {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN < 0.0)
        {
            sN = 0.0;
            tN = e;
            tD = c;//keeping the distance of closest approach perpendicular to the exiting track
                   //will later check if that's feasible consdiering the second segment length
        }
        else if (sN > sD)
        {
            sN = sD;
            tN = b+e;
            tD = c;//keeping the distance of closest approach perpendicular to the exiting track   
                   //will later check if that's feasible consdiering the second segment length
        }        
    }
    if (tN < 0.0)
    {
        tN = 0.0;
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    {
        tN = tD;
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else
        {
            sN = (-d + b);
            sD = a;
        }
    }
    if (fabs(sN) < pow(10.0,-8.0)) sN = 0.0;
    if (fabs(tN) < pow(10.0,-8.0)) tN = 0.0;
    dist = (point1 - point2) + (sN/sD)*v1 - (tN/tD)*v2;
    z = (point1 + (sN/sD)*v1)(2);
    return dist.Mag(); 
}

void SmearIncomingVectorAndEntrancePoint(Float_t& vec_x, Float_t& vec_y, Float_t& vec_z, Float_t& pos_x, Float_t& pos_y, Float_t& pos_z,
                                         Double_t thrown_x1, Double_t thrown_y1, Double_t thrown_x2, Double_t thrown_y2)
{
    //emulsion thickness in cm
    Double_t t = 0.02; //only plastic base
    TVector3 v(vec_x,vec_y,vec_z);
    TVector3 p(pos_x,pos_y,pos_z);
    TVector3 n(0.,0.,1.);
    TVector3 r0(0.,0.,-1.*t);
    TVector3 entering;
    //calculating the point where the hadron eneters the emulsion
    entering = p + (r0-p).Dot(n)/(v.Dot(n))*v;
    //smearing the enetering and exiting points (and hence the track direction)
    pos_x += thrown_x1;
    pos_y += thrown_y1;
    TVector3 exiting(pos_x,pos_y,pos_z);
    TVector3 smearing_entering(thrown_x2,thrown_y2,0.);
    entering += smearing_entering;
    //calculating the smeared track direction
    vec_x = exiting(0)-entering(0);
    vec_y = exiting(1)-entering(1);
}

void SmearOutgoingVectorAndExitingPoint(Float_t& vec_x, Float_t& vec_y, Float_t& vec_z, Float_t& pos_x, Float_t& pos_y, Float_t& pos_z,
                                         Double_t thrown_x1, Double_t thrown_y1, Double_t thrown_x2, Double_t thrown_y2)
{
    //emulsion thickness in cm
    Double_t t = 0.02; //only plastic base
    TVector3 v(vec_x,vec_y,vec_z);
    TVector3 p(pos_x,pos_y,pos_z);
    TVector3 n(0.,0.,1.);
    TVector3 r0(0.,0.,2.+t);
    TVector3 exiting;
    //calculating the point where the hadron exits the emulsion
    exiting = p + (r0-p).Dot(n)/(v.Dot(n))*v;
    //smearing the enetering and exiting points (and hence the track direction)
    pos_x += thrown_x1;
    pos_y += thrown_y1;
    TVector3 entering(pos_x,pos_y,pos_z);
    TVector3 smearing_exiting(thrown_x2,thrown_y2,0.);
    exiting += smearing_exiting;
    //calculating the smeared track direction
    vec_x = exiting(0)-entering(0);
    vec_y = exiting(1)-entering(1);
}

//this has now been tested!
void SmearVectorByDetectorResoultion(Float_t& vec_x, Float_t& vec_y, Float_t& vec_z, Double_t theta, Double_t phi)
{
    //rotate coordinate system until z-axis coincides with the vector that is to be smeared
    //use Euler angles
    //single vector describes 2 Euler angles (unless its magnitude is also interpreted as an angle)
    //alpha Euler angle
    Double_t alpha;
    if (vec_y>=0.) { alpha = acos(vec_x/sqrt(pow(vec_x,2.)+pow(vec_y,2.))); }
    else { alpha = 2.*PI - acos(vec_x/sqrt(pow(vec_x,2.)+pow(vec_y,2.))); }
    //rotation around z-axis by Euler angle alpha
    Double_t Rot_z_alpha_elements[9] = {cos(alpha),sin(alpha),0.,-1.*sin(alpha),cos(alpha),0.,0.,0.,1.};
    TMatrixD Rot_z_alpha(3,3,Rot_z_alpha_elements);
    //beta Euler angle
    Double_t beta = acos(vec_z);
    //rotation around y'-axis by Euler angle beta
    Double_t Rot_y_prime_beta_elements[9] = {cos(beta),0.,-1.*sin(beta),0.,1.,0.,sin(beta),0.,cos(beta)};
    TMatrixD R(3,3,Rot_y_prime_beta_elements);
    R *= Rot_z_alpha; //R is now the matrix for going from the initial to the rotated coordinate system
    //in the rotated coordinate system our unit vector points along the new z-axis
    //now we want to set its coordinates to (r=1,theta,phi) to account for the effect of smearing
    Double_t Unit_vector_elements[3] = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
    TMatrixD Unit_vector(3,1,Unit_vector_elements);
    //calculate the smeared momentum vector in our original coordiante system
    //first invert matrix R
    R.InvertFast();
    Unit_vector = R*Unit_vector;
    //set output values for the smeared vector coordinates
    vec_x = Unit_vector(0,0);
    vec_y = Unit_vector(1,0);
    vec_z = Unit_vector(2,0);
}

int main(int argc, char *argv[])
{
    char *thinmcdata = NULL;
    char *outfil = NULL;
    string thinarg("-thin");
    string foutarg("-o");
    ProcessDir dpr;
    vector<string> thinfiles;
    double files_thin = 0.0;
    for (int i=1; i<argc; i++) {
        if (thinarg.compare(argv[i])==0 && !thinmcdata) {
            thinmcdata = argv[++i];
            cout << "Thin target MC data is in  " << thinmcdata << endl;
            if (dpr.GetRootFiles(thinmcdata,thinfiles)==0) {
                cout << "Found " << thinfiles.size() << " file(s) in this directory" << endl;
                files_thin = (double)thinfiles.size();
            }
            else
            {
                cout << "Error occured";
                exit(-1);
            }
        }
        if (foutarg.compare(argv[i])==0 && !outfil) {
            outfil = argv[++i];
            cout << "Output file is set to " << outfil << endl;
        }
    }
    if (!outfil || !thinmcdata) {
        cerr << "Usage:" <<endl;
        cerr << "The input dir has to be specified" << endl;
        cerr << " -thin   <Thin MC data dir>" << endl;
        cerr << "The output file name has to be specified" << endl;
        cerr << " -o   <output file name>" << endl;
        exit(-1);
    }
    TFile fout(outfil,"RECREATE");
    Int_t nentries_thin;
    TChain *fChain;
    Int_t fCurrent;
    //Declaration of leaf types
    Int_t protnum;
    Int_t ipart;
    Float_t mom;
    Float_t pos[3];
    Float_t vec[3];
    Int_t pgen;
    Int_t ng;
    Float_t gpx[20]; 
    Float_t gpy[20]; 
    Float_t gpz[20]; 
    Float_t gvx[20]; 
    Float_t gvy[20]; 
    Float_t gvz[20]; 
    Int_t gpid[20]; 
    Int_t gmec[20]; 
    Bool_t flag;
    Float_t mom_enter;
    Float_t pos_enter[3];
    Float_t vec_enter[3];
    //List of branches
    TBranch *b_protnum;
    TBranch *b_ipart;
    TBranch *b_mom;
    TBranch *b_pos;
    TBranch *b_vec;
    TBranch *b_pgen;
    TBranch *b_ng;
    TBranch *b_gpx;
    TBranch *b_gpy;
    TBranch *b_gpz;
    TBranch *b_gvx;
    TBranch *b_gvy;
    TBranch *b_gvz;
    TBranch *b_gpid;
    TBranch *b_gmec;
    TBranch *b_flag;
    TBranch *b_mom_enter;
    TBranch *b_pos_enter;
    TBranch *b_vec_enter;
    //Initialise the TChain which will contain all thin target root files
    if (files_thin) {
        fChain = new TChain("h1000");
        fChain->Add(thinfiles[0].c_str());
        fCurrent = -1;
        fChain->SetMakeClass(1);
        fChain->SetBranchAddress("protnum", &protnum, &b_protnum);
        fChain->SetBranchAddress("ipart", &ipart, &b_ipart);
        fChain->SetBranchAddress("mom", &mom, &b_mom);
        fChain->SetBranchAddress("pos", &pos, &b_pos);
        fChain->SetBranchAddress("vec", &vec, &b_vec);
        fChain->SetBranchAddress("pgen", &pgen, &b_pgen);
        fChain->SetBranchAddress("ng", &ng, &b_ng);
        fChain->SetBranchAddress("gpx", &gpx, &b_gpx);
        fChain->SetBranchAddress("gpy", &gpy, &b_gpy);
        fChain->SetBranchAddress("gpz", &gpz, &b_gpz);
        fChain->SetBranchAddress("gvx", &gvx, &b_gvx);
        fChain->SetBranchAddress("gvy", &gvy, &b_gvy);
        fChain->SetBranchAddress("gvz", &gvz, &b_gvz);
        fChain->SetBranchAddress("gpid", &gpid, &b_gpid);
        fChain->SetBranchAddress("gmec", &gmec, &b_gmec);
        fChain->SetBranchAddress("flag", &flag, &b_flag);
        fChain->SetBranchAddress("mom_enter", &mom_enter, &b_mom_enter);
        fChain->SetBranchAddress("pos_enter", &pos_enter, &b_pos_enter);
        fChain->SetBranchAddress("vec_enter", &vec_enter, &b_vec_enter);
        for (int i=1; i<(int)files_thin; i++) {
            fChain->Add(thinfiles[i].c_str());
        }
        nentries_thin = fChain->GetEntries();
        cout << "Thin MC Entries " << nentries_thin << endl;
    }
    TH1D *h_matched_min_distance = new TH1D("matched_min_distance","matched_min_distance",140,0.,35.);
    TH1D *h_matched_z_coordinate = new TH1D("matched_z_coordinate","matched_z_coordinate",100,-1.,2.);
    TH1D *h_correct_matching_fraction = new TH1D("correct_matching_fraction","correct matching fraction",200,0.,0.5);
    TH1D *h_correct_matching_fraction_clone = (TH1D*)h_correct_matching_fraction->Clone("correct_matching_fraction_clone");
    TH1D *h_matched_min_distance_no_prot = new TH1D("matched_min_distance_no_prot","matched_min_distance_no_prot",140,0.,35.);
    TH1D *h_outgoing_theta = new TH1D("outgoing_theta","outgoing_theta",200,0.,PI);
    TH1D *h_truth_proton_matching_pattern = new TH1D("truth_proton_matching_pattern","truth_proton_matching_pattern",7,1.0,8.0);
    TH1D *h_theta_smearing = new TH1D("theta_smearing","theta_smearing",200,0.,2.);
    TH1D *h_phi_smearing = new TH1D("phi_smearing","phi_smearing",200,0.,PI/2.);
    TH2D *h_outgoing_proton = new TH2D("outgoing_proton","outgoing_proton",150,0.,5.,100,0.,2.);
    TH2D *h_smeared_outgoing_proton = new TH2D("smeared_outgoing_proton","smeared_outgoing_proton",125,0.,25.,50,0.,2.);
    //position where the proton is incident on the target upstream face
    vector<Float_t> vector_of_incoming_posx, vector_of_incoming_posy;
    //incoming momentum direction (unit vector) of the proton incident on the target
    vector<Float_t> vector_of_incoming_vecx, vector_of_incoming_vecy, vector_of_incoming_vecz;
    //incoming proton event id
    vector<Int_t> vector_of_incoming_event_id;
    //exiting position of hadrons leaving the target
    vector<Float_t> vector_of_outgoing_posx, vector_of_outgoing_posy, vector_of_outgoing_posz;
    //outgoing momentum direction (unit vector) of hadrons leaving the target
    vector<Float_t> vector_of_outgoing_vecx, vector_of_outgoing_vecy, vector_of_outgoing_vecz;
    //outgoing hadron id and the associated event id (same for all secondary hadrons originating from the same primary proton)
    vector<Int_t> vector_of_outgoing_pid, vector_of_outgoing_event_id;

    vector_of_incoming_posx.reserve(531000);
    vector_of_incoming_posy.reserve(531000);
    vector_of_incoming_vecx.reserve(531000);
    vector_of_incoming_vecy.reserve(531000);
    vector_of_incoming_vecz.reserve(531000);
    vector_of_incoming_event_id.reserve(531000);
    vector_of_outgoing_posx.reserve(nentries_thin);
    vector_of_outgoing_posy.reserve(nentries_thin);
    vector_of_outgoing_posz.reserve(nentries_thin);
    vector_of_outgoing_vecx.reserve(nentries_thin);
    vector_of_outgoing_vecy.reserve(nentries_thin);
    vector_of_outgoing_vecz.reserve(nentries_thin);
    vector_of_outgoing_pid.reserve(nentries_thin);
    vector_of_outgoing_event_id.reserve(nentries_thin);

    Double_t thrown_x1 = 0.;
    Double_t thrown_y1 = 0.;
    Double_t thrown_x2 = 0.;
    Double_t thrown_y2 = 0.;
    Double_t res = pow(5.,-5.);
    double dummy_dir[3];
    Int_t good_nentries_thin = 0;
    Int_t good_proton_nentries_thin = 0;
    TRandom2 rand;
    rand.SetSeed(1397);
    if (files_thin)
    {
        for (int i=0; i<nentries_thin; i++) {
            if (fChain) {
                Long64_t centry = fChain->LoadTree((Long64_t)i);
                if (!(centry<0) && fChain->InheritsFrom(TChain::Class())) {
                    TChain *chain = (TChain*)fChain;
                    if (chain->GetTreeNumber() != fCurrent) {
                        fCurrent = chain->GetTreeNumber();
                    }
                }
            }
            fChain->GetEntry((Long64_t)i);
            if (vec[2]<0. || vec_enter[2]<0.) continue; //bad MC entries cannot pass this criterion
            good_nentries_thin++;
            if (gpid[ng-1]==14) {
                good_proton_nentries_thin++;
                Double_t dummy_z = -1000.;
                Double_t dummy_closest_approach_distance = DistanceAtPointOfClosestApproach_LineSegments(vec_enter[0],vec_enter[1],vec_enter[2],
                                                                                            pos_enter[0],pos_enter[1],pos_enter[2],
                                                                                            vec[0],vec[1],vec[2],
                                                                                            pos[0],pos[1],pos[2],dummy_z);
                Double_t dummy_theta = acos(vec[2]/pow(pow(vec[0],2.)+pow(vec[1],2.)+pow(vec[2],2.),.5));
                h_outgoing_proton->Fill(dummy_closest_approach_distance*pow(10.,4.),dummy_theta*pow(10.,3.));
            }

            //Apply emulsion detector smearing on incoming and outgoing hadrons
            dummy_dir[0] = vec[0];
            dummy_dir[1] = vec[1];
            dummy_dir[2] = vec[2];
            //throw x and y positions to account for emulsion resoultion
            thrown_x1 = rand.Gaus(0.,res);
            thrown_y1 = rand.Gaus(0.,res);
            thrown_x2 = rand.Gaus(0.,res);
            thrown_y2 = rand.Gaus(0.,res);
            SmearOutgoingVectorAndExitingPoint(vec[0],vec[1],vec[2],
                                               pos[0],pos[1],pos[2],
                                               thrown_x1,thrown_y1,thrown_x2,thrown_y2);
            if (fabs(dummy_dir[0]*vec[0]+dummy_dir[1]*vec[1]+dummy_dir[2]*vec[2])<=1.) {
                h_theta_smearing->Fill(acos(dummy_dir[0]*vec[0]+dummy_dir[1]*vec[1]+dummy_dir[2]*vec[2])*1000);
            }
            if (fabs(dummy_dir[0]*vec[0]+dummy_dir[1]*vec[1])/(pow(pow(dummy_dir[0],2.)+pow(dummy_dir[1],2.),.5)*pow(pow(vec[0],2.)+pow(vec[1],2.),.5))) {
                h_phi_smearing->Fill(acos((dummy_dir[0]*vec[0]+dummy_dir[1]*vec[1])/(pow(pow(dummy_dir[0],2.)+pow(dummy_dir[1],2.),.5)*pow(pow(vec[0],2.)+pow(vec[1],2.),.5))));
            } 
            
            vector_of_outgoing_posx.push_back(pos[0]);
            vector_of_outgoing_posy.push_back(pos[1]);
            vector_of_outgoing_posz.push_back(pos[2]);
            vector_of_outgoing_vecx.push_back(vec[0]);
            vector_of_outgoing_vecy.push_back(vec[1]);
            vector_of_outgoing_vecz.push_back(vec[2]);
            vector_of_outgoing_pid.push_back(gpid[ng-1]);
            vector_of_outgoing_event_id.push_back(protnum);

            //throw x and y positions to account for emulsion resoultion
            thrown_x1 = rand.Gaus(0.,res);
            thrown_y1 = rand.Gaus(0.,res);
            thrown_x2 = rand.Gaus(0.,res);
            thrown_y2 = rand.Gaus(0.,res);
            SmearIncomingVectorAndEntrancePoint(vec_enter[0],vec_enter[1],vec_enter[2],
                                                pos_enter[0],pos_enter[1],pos_enter[2],
                                                thrown_x1,thrown_y1,thrown_x2,thrown_y2);
            //cout << "Event ID Associated with exiting hadron:" << protnum << endl;
            //check if the smeared outgoing proton histogram needs to be filled
            if (gpid[ng-1]==14) {
                Double_t dummy_z = -1000.;
                Double_t dummy_closest_approach_distance = DistanceAtPointOfClosestApproach_LineSegments(vec_enter[0],vec_enter[1],vec_enter[2],
                                                                                            pos_enter[0],pos_enter[1],pos_enter[2],
                                                                                            vec[0],vec[1],vec[2],
                                                                                            pos[0],pos[1],pos[2],dummy_z);
                Double_t dummy_theta = acos(vec[2]/pow(pow(vec[0],2.)+pow(vec[1],2.)+pow(vec[2],2.),.5));
                h_smeared_outgoing_proton->Fill(dummy_closest_approach_distance*pow(10.,4.),dummy_theta*pow(10.,3.)); 
            }
            if (i>0 && vector_of_incoming_event_id.back()==protnum) continue;
            dummy_dir[0] = vec_enter[0];
            dummy_dir[1] = vec_enter[1];
            dummy_dir[2] = vec_enter[2];
            if (fabs(dummy_dir[0]*vec_enter[0]+dummy_dir[1]*vec_enter[1]+dummy_dir[2]*vec_enter[2])>=1.) {
                h_theta_smearing->Fill(acos(dummy_dir[0]*vec_enter[0]+dummy_dir[1]*vec_enter[1]+dummy_dir[2]*vec_enter[2])*1000);
            }
            if (fabs(dummy_dir[0]*vec_enter[0]+dummy_dir[1]*vec_enter[1])/(pow(pow(dummy_dir[0],2.)+pow(dummy_dir[1],2.),.5)*pow(pow(vec_enter[0],2.)+pow(vec_enter[1],2.),.5))) {
                h_phi_smearing->Fill(acos((dummy_dir[0]*vec_enter[0]+dummy_dir[1]*vec_enter[1])/(pow(pow(dummy_dir[0],2.)+pow(dummy_dir[1],2.),.5)*pow(pow(vec_enter[0],2.)+pow(vec_enter[1],2.),.5))));
            } 
            vector_of_incoming_posx.push_back(pos_enter[0]);
            vector_of_incoming_posy.push_back(pos_enter[1]);
            vector_of_incoming_vecx.push_back(vec_enter[0]);
            vector_of_incoming_vecy.push_back(vec_enter[1]);
            vector_of_incoming_vecz.push_back(vec_enter[2]);
            vector_of_incoming_event_id.push_back(protnum);

            //cout << "Primary Proton Event ID:" << protnum << endl;
        }
    }
    cout << "Size of primary proton Event ID array:" << vector_of_incoming_event_id.size() << endl;
    cout << "Size of exiting hadron Event ID array:" << vector_of_outgoing_event_id.size() << endl;
    cout << "Starting to match incoming and outgoing tracks" << endl;
    Int_t good_matches = 0;
    Int_t good_proton_matches = 0;
    //every exiting particle direction can be matched to an incoming track based on the minimum distance
    std::vector<Double_t> min_distance;
    min_distance.reserve(nentries_thin);
    std::vector<Double_t> min_distance_2nd;
    min_distance_2nd.reserve(nentries_thin);
    std::vector<Double_t> min_distance_3rd;
    min_distance_3rd.reserve(nentries_thin);
    std::vector<Double_t> min_distance_4th;
    min_distance_4th.reserve(nentries_thin);
    std::vector<Double_t> min_distance_5th;
    min_distance_5th.reserve(nentries_thin);
    std::vector<Double_t> min_distance_6th;
    min_distance_6th.reserve(nentries_thin);
    std::vector<Double_t> min_distance_7th;
    min_distance_7th.reserve(nentries_thin);

    std::vector<Int_t> matched_incoming_track;
    matched_incoming_track.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_2nd;
    matched_incoming_track_2nd.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_3rd;
    matched_incoming_track_3rd.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_4th;
    matched_incoming_track_4th.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_5th;
    matched_incoming_track_5th.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_6th;
    matched_incoming_track_6th.reserve(nentries_thin);
    std::vector<Int_t> matched_incoming_track_7th;
    matched_incoming_track_7th.reserve(nentries_thin); 

    std::vector<Double_t> matched_z_coordinate;
    matched_z_coordinate.reserve(nentries_thin);
    std::vector<Int_t> correct_matching_trigger;
    correct_matching_trigger.reserve(nentries_thin);
    //time_t start, end;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    if (files_thin)
    {
        for (int i=0; i<good_nentries_thin; i++) {
            if ((i%1000)==0) { cout << "Now matching the exiting track for the " << i << "th hadron..."  << endl; }
            //cout << "Now matching the exiting track for the " << i << "th hadron..."  << endl;
            min_distance.push_back(pow(10.,5.));
            min_distance_2nd.push_back(pow(10.,5.));
            min_distance_3rd.push_back(pow(10.,5.));
            min_distance_4th.push_back(pow(10.,5.));
            min_distance_5th.push_back(pow(10.,5.));
            min_distance_6th.push_back(pow(10.,5.));
            min_distance_7th.push_back(pow(10.,5.));

            matched_incoming_track.push_back(-1);
            matched_incoming_track_2nd.push_back(-1);
            matched_incoming_track_3rd.push_back(-1);
            matched_incoming_track_4th.push_back(-1);
            matched_incoming_track_5th.push_back(-1);
            matched_incoming_track_6th.push_back(-1);
            matched_incoming_track_7th.push_back(-1);

            matched_z_coordinate.push_back(-100.);
            correct_matching_trigger.push_back(0);
            Float_t vec2_x = vector_of_outgoing_vecx[i];
            Float_t vec2_y = vector_of_outgoing_vecy[i];
            Float_t vec2_z = vector_of_outgoing_vecz[i];
            Float_t p2_x = vector_of_outgoing_posx[i];
            Float_t p2_y = vector_of_outgoing_posy[i];
            Float_t p2_z = vector_of_outgoing_posz[i];
            for (int j=0; j<vector_of_incoming_event_id.size(); j++) {
                Float_t vec1_x = vector_of_incoming_vecx[j];
                Float_t vec1_y = vector_of_incoming_vecy[j];
                Float_t vec1_z = vector_of_incoming_vecz[j];
                Float_t p1_x = vector_of_incoming_posx[j];
                Float_t p1_y = vector_of_incoming_posy[j];
                Float_t p1_z = 0.;
                if (i==0 && j==0) 
                { 
                    //time(&start);
                    start = std::chrono::system_clock::now();
                    std::time_t now_c = std::chrono::system_clock::to_time_t(start); 
                    cout << "DistanceAtPointOfClosestApproach_LineSegments started running at: " << ctime(&now_c);
                }
                Double_t potential_matched_z_coordinate = -1000.;
                Double_t potential_min_distance = DistanceAtPointOfClosestApproach_LineSegments(vec1_x, vec1_y, vec1_z, p1_x, p1_y, p1_z, vec2_x, vec2_y, vec2_z, p2_x, p2_y, p2_z, potential_matched_z_coordinate); 
                if (i==0 && j==0) 
                { 
                    //time(&end);
                    end = std::chrono::system_clock::now();
                    std::chrono::duration<double> elapsed_seconds = end-start;
                    cout << "Elapsed time while running DistanceAtPointOfClosestApproach_LineSegments is:" << elapsed_seconds.count() << endl;
                    //cout << "DistanceAtPointOfClosestApproach_LineSegments finished running at: " << ctime(&end);
                }
                
                if (potential_min_distance < min_distance[i])
                {
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = min_distance_5th[i];
                    matched_incoming_track_6th[i] = matched_incoming_track_5th[i];
                    min_distance_5th[i] = min_distance_4th[i];
                    matched_incoming_track_5th[i] = matched_incoming_track_4th[i];
                    min_distance_4th[i] = min_distance_3rd[i];
                    matched_incoming_track_4th[i] = matched_incoming_track_3rd[i];
                    min_distance_3rd[i] = min_distance_2nd[i];
                    matched_incoming_track_3rd[i] = matched_incoming_track_2nd[i];
                    min_distance_2nd[i] = min_distance[i];
                    matched_incoming_track_2nd[i] = matched_incoming_track[i];
                    min_distance[i] = potential_min_distance;
                    matched_incoming_track[i] = j;
                    
                    matched_z_coordinate[i] = potential_matched_z_coordinate;
                }
                else if (potential_min_distance < min_distance_2nd[i])
                {
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = min_distance_5th[i];
                    matched_incoming_track_6th[i] = matched_incoming_track_5th[i];
                    min_distance_5th[i] = min_distance_4th[i];
                    matched_incoming_track_5th[i] = matched_incoming_track_4th[i];
                    min_distance_4th[i] = min_distance_3rd[i];
                    matched_incoming_track_4th[i] = matched_incoming_track_3rd[i];
                    min_distance_3rd[i] = min_distance_2nd[i];
                    matched_incoming_track_3rd[i] = matched_incoming_track_2nd[i];
                    min_distance_2nd[i] = potential_min_distance;
                    matched_incoming_track_2nd[i] = j;
                }
                else if (potential_min_distance < min_distance_3rd[i])
                {
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = min_distance_5th[i];
                    matched_incoming_track_6th[i] = matched_incoming_track_5th[i];
                    min_distance_5th[i] = min_distance_4th[i];
                    matched_incoming_track_5th[i] = matched_incoming_track_4th[i];
                    min_distance_4th[i] = min_distance_3rd[i];
                    matched_incoming_track_4th[i] = matched_incoming_track_3rd[i];
                    min_distance_3rd[i] = potential_min_distance;
                    matched_incoming_track_3rd[i] = j;
                }
                else if (potential_min_distance < min_distance_4th[i])
                {   
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = min_distance_5th[i];
                    matched_incoming_track_6th[i] = matched_incoming_track_5th[i];
                    min_distance_5th[i] = min_distance_4th[i];
                    matched_incoming_track_5th[i] = matched_incoming_track_4th[i];
                    min_distance_4th[i] = potential_min_distance;
                    matched_incoming_track_4th[i] = j;
                }
                else if (potential_min_distance < min_distance_5th[i])
                {   
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = min_distance_5th[i];
                    matched_incoming_track_6th[i] = matched_incoming_track_5th[i];
                    min_distance_5th[i] = potential_min_distance;
                    matched_incoming_track_5th[i] = j;
                }
                else if (potential_min_distance < min_distance_6th[i])
                {
                    min_distance_7th[i] = min_distance_6th[i];
                    matched_incoming_track_7th[i] = matched_incoming_track_6th[i];
                    min_distance_6th[i] = potential_min_distance;
                    matched_incoming_track_6th[i] = j;
                }
                else if (potential_min_distance < min_distance_7th[i])
                {
                    min_distance_7th[i] = potential_min_distance;
                    matched_incoming_track_7th[i] = j;
                }
            }
            //cout << "Minimum distance between the exiting hadron's track and one of incoming tracks is:" << min_distance[i] << endl;
            //cout << "Matched to " << matched_incoming_track[i] << "th incoming proton"  << endl;
            if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track[i])]) 
            { 
                good_matches++;
                if (vector_of_outgoing_pid[i] == 14) {
                    good_proton_matches++;
                }
                if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(1.0);
                correct_matching_trigger[i] = 1;
                //cout << "SUCCESS in matching!" << endl; 
            }
            else
            {
                if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_2nd[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(2.0);
                }
                else if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_3rd[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(3.0);
                }
                else if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_4th[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(4.0);
                }
                else if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_5th[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(5.0);
                }
                else if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_6th[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(6.0);
                }
                else if (vector_of_outgoing_event_id[i] == vector_of_incoming_event_id[(matched_incoming_track_7th[i])])
                {
                    if (vector_of_outgoing_pid[i]==14) h_truth_proton_matching_pattern->Fill(7.0);
                }
            }
            h_outgoing_theta->Fill(acos(vector_of_outgoing_vecz[i]));
            h_matched_min_distance->Fill(min_distance[i]*pow(10.,4.));
            h_matched_z_coordinate->Fill(matched_z_coordinate[i]);
            Double_t dummy_theta = acos(vector_of_outgoing_vecz[i]/pow(pow(vector_of_outgoing_vecx[i],2.)+pow(vector_of_outgoing_vecy[i],2.)+pow(vector_of_outgoing_vecz[i],2.),.5));
            if (dummy_theta<0.5) {
                h_correct_matching_fraction->Fill(dummy_theta,correct_matching_trigger[i]);
                h_correct_matching_fraction_clone->Fill(dummy_theta);
            }
            if (vector_of_outgoing_pid[i]!=14) { h_matched_min_distance_no_prot->Fill(min_distance[i]*pow(10.,4.)); }
        }
    }
    for (int i=1; i<=h_correct_matching_fraction->GetNbinsX(); i++) {
        if (h_correct_matching_fraction_clone->GetBinContent(i)==0) continue;
        Double_t dummy_fraction = h_correct_matching_fraction->GetBinContent(i)/h_correct_matching_fraction_clone->GetBinContent(i);
        h_correct_matching_fraction->SetBinContent(i,dummy_fraction);
    }
    h_correct_matching_fraction->GetXaxis()->SetTitle("theta[rad]");
    h_correct_matching_fraction->GetYaxis()->SetTitle("Fraction");
    h_correct_matching_fraction->Write();
    h_correct_matching_fraction_clone->Write();
    h_theta_smearing->GetXaxis()->SetTitle("theta[mrad]");
    h_theta_smearing->Write();
    h_phi_smearing->GetXaxis()->SetTitle("phi[rad]");
    h_phi_smearing->Write();    
    h_matched_min_distance->Write();
    h_matched_z_coordinate->Write();
    h_matched_min_distance_no_prot->Write();
    h_outgoing_theta->GetXaxis()->SetTitle("theta[rad]");
    h_outgoing_theta->Write();
    h_truth_proton_matching_pattern->Write();
    h_outgoing_proton->GetXaxis()->SetTitle("distance between tracks[micron]");
    h_outgoing_proton->GetYaxis()->SetTitle("exiting theta[mrad]");
    h_outgoing_proton->Write();
    h_smeared_outgoing_proton->GetXaxis()->SetTitle("distance between tracks[micron]");
    h_smeared_outgoing_proton->GetYaxis()->SetTitle("exiting theta[mrad]");
    h_smeared_outgoing_proton->Write();
    cout << "Total number of forward exiting hadrons:" << good_nentries_thin << endl;
    cout << "Fraction of exiting hadrons correctly matched to incident protons:" << good_matches*1./(good_nentries_thin*1.) << endl;
    cout << "Fraction of exiting protons correctly matched to incident protons:" << good_proton_matches*1./(good_proton_nentries_thin*1.) << endl;
    cout << "Number of forward exiting protons:" << good_proton_nentries_thin << endl;
    
    Int_t number_of_incoming_line_segments = vector_of_incoming_event_id.size();    
    TPolyLine3D *l_incoming[number_of_incoming_line_segments];
    TCanvas *c1 = new TCanvas("incoming_line_segments","incoming_line_segments",500,500);
    TView *view1 = TView::CreateView(1);
    for (int i=0; i<number_of_incoming_line_segments; i++) {
        l_incoming[i] = new TPolyLine3D(2);
        PrimaryProton(vector_of_incoming_vecx[i],vector_of_incoming_vecy[i],vector_of_incoming_vecz[i],vector_of_incoming_posx[i],vector_of_incoming_posy[i],0.);
        l_incoming[i]->SetPoint(0,vector_of_incoming_posx[i]-0.16,vector_of_incoming_posy[i]-0.21,0.);
        cout << "Line segment start:(" << vector_of_incoming_posx[i] << "," << vector_of_incoming_posy[i] << ",0.)" << endl;
        l_incoming[i]->SetPoint(1,vector_of_incoming_vecx[i]+vector_of_incoming_posx[i]-0.16,
                         vector_of_incoming_vecy[i]+vector_of_incoming_posy[i]-0.21,
                         vector_of_incoming_vecz[i]);
        cout << "Line segment end:(" << vector_of_incoming_vecx[i]+vector_of_incoming_posx[i]-0.16 
                              << "," << vector_of_incoming_vecy[i]+vector_of_incoming_posy[i]-0.21 
                              << "," << vector_of_incoming_vecz[i] << ")" << endl;
        l_incoming[i]->Draw();
    }
    view1->ShowAxis();
    c1->Write();

    Int_t number_of_outgoing_line_segments = good_nentries_thin;
    TPolyLine3D *l_outgoing[number_of_outgoing_line_segments];
    TCanvas *c2 = new TCanvas("outgoing_line_segments","outgoing_line_segments",500,500);
    TView *view2 = TView::CreateView(1);
    for (int i=0; i<number_of_outgoing_line_segments; i++) {
        if (vector_of_outgoing_pid[i]!=14) continue;
        l_outgoing[i] = new TPolyLine3D(2);
        ExitingHadron(vector_of_outgoing_vecx[i],vector_of_outgoing_vecy[i],vector_of_outgoing_vecz[i],vector_of_outgoing_posx[i],vector_of_outgoing_posy[i],vector_of_outgoing_posz[i]);
        l_outgoing[i]->SetPoint(0,vector_of_outgoing_posx[i]-0.16,vector_of_outgoing_posy[i]-0.21,vector_of_outgoing_posz[i]);
        cout << "Line segment start:(" << vector_of_outgoing_posx[i] << "," << vector_of_outgoing_posy[i] << "," << vector_of_outgoing_posz[i] << ")" << endl;
        if (pow(pow(vector_of_outgoing_posx[i]-0.16,2.)+pow(vector_of_outgoing_posy[i]-0.21,2.),.5)>target_radius || vector_of_outgoing_posz[i]>=0.) { cout << "Check point above!" << endl; }
        l_outgoing[i]->SetPoint(1,vector_of_outgoing_vecx[i]+vector_of_outgoing_posx[i]-0.16,
                                  vector_of_outgoing_vecy[i]+vector_of_outgoing_posy[i]-0.21,
                                  vector_of_outgoing_vecz[i]+vector_of_outgoing_posz[i]);
        cout << "Line segment end:(" << vector_of_outgoing_vecx[i]+vector_of_outgoing_posx[i]-0.16
                              << "," << vector_of_outgoing_vecy[i]+vector_of_outgoing_posy[i]-0.21
                              << "," << vector_of_outgoing_vecz[i]+vector_of_outgoing_posz[i] << ")" << endl;
        if (pow(pow(vector_of_outgoing_vecx[i]+vector_of_outgoing_posx[i]-0.16,2.)+pow(vector_of_outgoing_vecy[i]+vector_of_outgoing_posy[i]-0.21,2.),.5)>target_radius
            || vector_of_outgoing_vecz[i]<0.) { cout << "Check point above!" << endl; }
        l_outgoing[i]->Draw();
    }
    view2->ShowAxis();
    c2->Write();
    cout << "Done" <<endl;
    fout.Close();
    return 0;
}

