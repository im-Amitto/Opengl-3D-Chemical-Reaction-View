//Amit Meena
//160001004
//Pushpendra Kumar
//16001046
#include <bits/stdc++.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include "lib/RgbImage.cpp"
#include "lib/RgbImage.h"

using namespace std;

//Menu Parameters
static int window;
static int menu_id;
static int submenu_id;
static int value = 0; 

//Flags deciding what part to render and whay to leave
int start =0;
int textMode = 1;
int floorMode = 1;
int presentation = 0;
float atomSize = 0.4;
bool bondMode = true;
bool atomMode = true;

//Menu Function
void menu(int num){
  if(num == 0){
    glutDestroyWindow(window);
    exit(0);
  }else{
    value = num;
  }
  glutPostRedisplay();
} 

//Menu Creator
void createMenu(void){     
    submenu_id = glutCreateMenu(menu);
    glutAddMenuEntry("Increase", 40);
    glutAddMenuEntry("Decrease", 41);   
    menu_id = glutCreateMenu(menu);
    if(atomMode) //Atom Toggle
    {
        glutAddMenuEntry("Atom - Off",20);
    }else
    {
        glutAddMenuEntry("Atom - On",21);
    }
    if(bondMode) //Bond Toggle
    {
        glutAddMenuEntry("Bonds - Off",30);
    }else
    {
        glutAddMenuEntry("Bonds - On",31);
    }
    if(presentation != -1) //Presentation Toggle
    {
        glutAddMenuEntry("Presentation - Off",50);
    }else
    {
        glutAddMenuEntry("Presentation - On",51);
    }
    if(textMode) //Text Toggle
    {
        glutAddMenuEntry("Info - Off",60);
    }else
    {
        glutAddMenuEntry("Info - On",61);
    }
    if(floorMode) //Floor Toggle
    {
        glutAddMenuEntry("Floor - Off",70);
    }else
    {
        glutAddMenuEntry("Floor - On",71);
    }
    glutAddSubMenu("Atom Size", submenu_id);
    glutAddMenuEntry("Quit", 0);     
    glutAttachMenu(GLUT_RIGHT_BUTTON);
} 
//Symbol to Atom Number
enum Symbol {INVALID = -1,H = 1, He,Li, Be, B, C, N, O, F,Ne,Na, Mg, Al, Si, P, S, Cl, Ar,K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,Cs, Ba, Hf=72, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,Fr, Ra};

GLuint texture[120];

string atomInfoS = "";
string Instruction = "Instructions :-\nClick on atom to get it's info\nN --> Next Reaction\nP --> Prev Reaction\nW,A,S,D,Q,E --> Camera Movement on X,Y,Z\nZ + Mouse Move --> Rotate Window";
string projectInfo = "Chemical Reaction 3D View\nSubmitted by :- \n  Amit Kumar Meena - 160001004 \n  Pushpendra Kumar - 160001046\nSupervised by :- \n  Dr. Somnath Dey";

string AtomInfo[120];

//Initialize Atomic Info
void detail()
{
AtomInfo[1]="Hydrogen: \n1. It is a colourless,odourless and tasteless gas.\n2. It is the lightest gas known.\n3. It is only very slightly soluble in water.\n4. It can be liquefied under high pressure and at low temperature.";
AtomInfo[6] ="Carbon:\n1. It is a non-metallic element.\n2. It occurs both in free as well as combined state.\n3. Air also contain carbon as carbon-di-oxide.\n4. In free state it occurs as diamond,coal and graphite.\n5. Carbon forms hydrites known as hydrocarbon.";
AtomInfo[7] ="Nitrogen:\n1. It is a typical non-metal.\n2. It exists as diatomic molecule.\n3. It is highly electronegative element.\n4. The oxidation state of nitrogen varies from -3 to +5.\n5. Molecular nitrogen is called dinitrogen.";
AtomInfo[8] ="Oxygen:\n1. It is a non-metal.\n2. It is paramagnatic in nature.\n3. It is most abundant element in earth's crust.\n4. It is di-atomic in nature.\n5. Molecular oxygen is also called dioxygen.";
AtomInfo[16] = "Sulphur:\n1.)It is present in small proportion\n2.)It is used as disinfectant for destroying bacteria,fungi\n3.)It is used in vulcanisation of rubber\n4.)It is a constituent of medicines for skin diseases\n5.)It is used in manufacture of matches,fire-works,etc";
AtomInfo[17] = "Chlorine:\n1.)It is a non-metallic element\n2.)It is used in bleaching textile,yarn,paper,pulp\n3.)It is used in the sterilization of drinking water\n4.)It is used in the manufacture of vinyl chloride\n5.)It is used in preparing insecticides such as D.D.T. & B.H.C";
AtomInfo[35] =  "Bromine:\n1.)It's atomic radius is 101.4pm\n2.)It's ionisation energy is 1142 KJ/mol\n3.)It is reddish brown in colour\n4.)It makes bromo compounds in organic chemistry\n5.)It finds use in medicine\n6.)It is a liquid with obnoxious smel";
}

//Reaction Details Structure
struct rxn{
    string a; //Molecule 1
    string b; //Molecule 2
    string c; //Result Molecule
    string name; // Reaction name
};

rxn mole[5];
int moleI=0;

//Initializing Reaction
void getRxns(){
mole[0] = {"./XYZ/benzene.xyz","./XYZ/h2so4.xyz","./XYZ/benzenesulfonic_acid.xyz","Sulphonication Reaction"};
mole[1] = {"./XYZ/benzene.xyz","./XYZ/hydrogen.xyz","./XYZ/cyclohexane.xyz","Hydrogenation Reaction"};
mole[2] = {"./XYZ/benzene.xyz","./XYZ/bromine.xyz","./XYZ/bromobenzene.xyz","Bromination Reaction"};
mole[3] = {"./XYZ/c4h8.xyz","./XYZ/chcl3.xyz","./XYZ/c4h6.xyz","Simon-Craft Reaction"};
mole[4] = {"./XYZ/benzene.xyz","./XYZ/ch3cl.xyz","./XYZ/toulene.xyz","Friedel-Craft Alkylation"};
}

//Floor Vertices
float fVert[4][3] = {
    {-50.0,6.0, -50.0},
    {+50.0,6.0, -50.0},
    {+50.0,6.0, +50.0},
    {-50.0,6.0, +50.0}
};

//Structure for atom details
struct AtomDetail
{
    AtomDetail(Symbol symbol, string name, double radius, double bondRadius)
    : symbol(symbol), name(name), radius(radius), bondRadius(bondRadius)
    {}
    Symbol symbol; //Atom Symbol
    string name; //Atom Name
    double radius,bondRadius; //Atomic Radius && Atom Bond Radius
};


//Structure for atom
struct Atom
{
    Atom(Symbol type) : symbol(type) {}
    Symbol symbol; //Atomic Symbol
    double x,y,z; //Co-ordinates of atom
    double DistanceSquared(Atom other) //Distance between the current atom and the atom which is passed 
    {
        return (x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z);
    }
};

//Structure for molecule
struct Molecule 
{
    vector<Atom> atoms;
    vector<tuple<int, int>> bonds;
};

//Loading and returning texture for the filname
GLuint loadTextureFromFile(char *filename)
{   
    GLuint temp = 0;
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel(GL_FLAT);
   glEnable(GL_DEPTH_TEST);

   RgbImage theTexMap( filename );

    glGenTextures(1, &temp);		
	glBindTexture(GL_TEXTURE_2D, temp);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	
	glTexImage2D(GL_TEXTURE_2D, 0, 3, theTexMap.GetNumCols(), theTexMap.GetNumRows(), 0, GL_RGB, GL_UNSIGNED_BYTE, theTexMap.ImageData() );
    return temp;
}

//Initialize a periodic table
struct PeriodicTable
{
    //initialize
    PeriodicTable()
    {
      table.push_back(AtomDetail(Symbol::H, "H", 0.6, .37));
        table.push_back(AtomDetail(Symbol::He, "He", 1.9, .32));
        table.push_back(AtomDetail(Symbol::Li, "Li", 1.55, 1.34));
        table.push_back(AtomDetail(Symbol::Be, "Be", 1.12, .90));
        table.push_back(AtomDetail(Symbol::B, "B", 0.98, .82));
        table.push_back(AtomDetail(Symbol::C, "C", 0.91, .77));
        table.push_back(AtomDetail(Symbol::N, "N", 0.92, .75));
        table.push_back(AtomDetail(Symbol::O, "O", 0.80, .73));
        table.push_back(AtomDetail(Symbol::F, "F", 0.57, .71));
        table.push_back(AtomDetail(Symbol::Ne, "Ne", 0.51, .69));
        table.push_back(AtomDetail(Symbol::Na, "Na", 1.9, 1.54));
        table.push_back(AtomDetail(Symbol::Mg, "Mg", 1.6, 1.30));
        table.push_back(AtomDetail(Symbol::Al, "Al", 1.4, 1.118));
        table.push_back(AtomDetail(Symbol::Si, "Si", 1.32, 1.11));
        table.push_back(AtomDetail(Symbol::P, "P", 1.28, 1.06));
        table.push_back(AtomDetail(Symbol::S, "S", 1.6, 1.02));
        table.push_back(AtomDetail(Symbol::Cl, "Cl", 1.55, .99));
        table.push_back(AtomDetail(Symbol::Ar, "Ar", 0.88, .97));
        table.push_back(AtomDetail(Symbol::K, "K", 0.88, 1.96));
        table.push_back(AtomDetail(Symbol::Ca, "Ca", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::Sc, "Sc", 0.88, 1.44));
        table.push_back(AtomDetail(Symbol::Ti, "Ti", 0.88, 1.36));
        table.push_back(AtomDetail(Symbol::V, "V", 0.88, 1.25));
        table.push_back(AtomDetail(Symbol::Cr, "Cr", 0.88, 1.27));
        table.push_back(AtomDetail(Symbol::Mn, "Mn", 0.88, 1.39));
        table.push_back(AtomDetail(Symbol::Fe, "Fe", 0.88, 1.25));
        table.push_back(AtomDetail(Symbol::Co, "Co", 0.88, 1.26));
        table.push_back(AtomDetail(Symbol::Ni, "Ni", 0.88, 1.21));
        table.push_back(AtomDetail(Symbol::Cu, "Cu", 0.88, 1.38));
        table.push_back(AtomDetail(Symbol::Zn, "Zn", 0.88, 1.31));
        table.push_back(AtomDetail(Symbol::Ga, "Ga", 0.88, 1.26));
        table.push_back(AtomDetail(Symbol::Ge, "Ge", 0.88, 1.22));
        table.push_back(AtomDetail(Symbol::As, "As", 0.88, 1.19));
        table.push_back(AtomDetail(Symbol::Se, "Se", 0.88, 1.16));
        table.push_back(AtomDetail(Symbol::Br, "Br", 1.45, 1.14));
        table.push_back(AtomDetail(Symbol::Kr, "Kr", 0.88, 1.10));
        table.push_back(AtomDetail(Symbol::Rb, "Rb", 0.88, 2.11));
        table.push_back(AtomDetail(Symbol::Sr, "Sr", 0.88, 1.92));
        table.push_back(AtomDetail(Symbol::Y, "Y", 0.88, 1.62));
        table.push_back(AtomDetail(Symbol::Zr, "Zr", 0.88, 1.48));
        table.push_back(AtomDetail(Symbol::Nb, "Nb", 0.88, 1.37));
        table.push_back(AtomDetail(Symbol::Mo, "Mo", 0.88, 1.45));
        table.push_back(AtomDetail(Symbol::Tc, "Tc", 0.88, 1.56));
        table.push_back(AtomDetail(Symbol::Ru, "Ru", 0.88, 1.26));
        table.push_back(AtomDetail(Symbol::Rh, "Rh", 0.88, 1.35));
        table.push_back(AtomDetail(Symbol::Pd, "Pd", 0.88, 1.31));
        table.push_back(AtomDetail(Symbol::Ag, "Ag", 0.88, 1.53));
        table.push_back(AtomDetail(Symbol::Cd, "Cd", 0.88, 1.48));
        table.push_back(AtomDetail(Symbol::In, "In", 0.88, 1.44));
        table.push_back(AtomDetail(Symbol::Sn, "Sn", 0.88, 1.41));
        table.push_back(AtomDetail(Symbol::Sb, "Sb", 0.88, 1.38));
        table.push_back(AtomDetail(Symbol::Te, "Te", 0.88, 1.35));
        table.push_back(AtomDetail(Symbol::I, "I", 0.88, 1.33));
        table.push_back(AtomDetail(Symbol::Xe, "Xe", 0.88, 1.30));
        table.push_back(AtomDetail(Symbol::Cs, "Cs", 0.88, 2.25));
        table.push_back(AtomDetail(Symbol::Ba, "Ba", 0.88, 1.98));
        table.push_back(AtomDetail(Symbol::Hf, "Hf", 0.88, 1.69));
        table.push_back(AtomDetail(Symbol::Ta, "Ta", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::W, "Re", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Os, "Os", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Ir, "Ir", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Pt, "Ru", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Au, "Au", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Hg, "Hg", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Ti, "Ti", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Pb, "Pb", 0.88, 1.0));
        table.push_back(AtomDetail(Symbol::Bi, "Bi", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::Po, "Po", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::At, "At", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::Rn, "Rn", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::Fr, "Fr", 0.88, 1.74));
        table.push_back(AtomDetail(Symbol::Ra, "Ra", 0.88, 1.74));
    }
    //check if two atom can bond with each other
    bool AtomsCanBond(Atom atom1, Atom atom2)
    {
        double distance = atom1.DistanceSquared(atom2);
        auto it1 = find_if(table.cbegin(), table.cend(), [&atom1](AtomDetail x) {
            return x.symbol == atom1.symbol;
        }); 
        auto it2 = find_if(table.cbegin(), table.cend(), [&atom2](AtomDetail x) {
            return x.symbol == atom2.symbol;
        }); 
        if (it1 == table.cend() || it2 == table.cend())
            return false;  
        double radius1 = (*it1).bondRadius;
        double radius2 = (*it2).bondRadius;
        double radius =  (radius1+radius2)*(radius1+radius2) * 1.1;
        if (distance <= radius) 
            return true;
        return false;
    }
    vector<AtomDetail> table;
};

//parse the file data to extract info from XYZ format and to store them in molecule details which is being used to render atom
struct DataParser {
    DataParser(PeriodicTable table)
    : chem(table) {}    
    Molecule ParseData(istream& stream)
    {
        Molecule molecule;
        string line;
        char atomName[64];
        
        while (getline(stream, line))
        {
            int count;
            double x, y, z;
            
            if (sscanf(line.c_str(), "%s %lf %lf %lf", atomName, &x, &y, &z) == 4)
            {
                auto it = find_if(chem.table.cbegin(), chem.table.cend(), [&atomName](AtomDetail x) {
                    return x.name == atomName;
                });
 
                if (it == chem.table.cend())
                    continue; 
                
                Atom newAtom((*it).symbol);
                newAtom.x = x;
                newAtom.y = y;
                newAtom.z = z;
                
                molecule.atoms.push_back(newAtom);
            }
            else if (sscanf(line.c_str(), "%i", &count) == 1)
            {
                if (count > 0)
                {
                    molecule.atoms.reserve(count);
                }
            }
        }
        
        
        for (auto it = molecule.atoms.begin(); it != molecule.atoms.end(); it ++)
        {
            Atom atom = (*it);
            
            for (auto j = molecule.atoms.begin(); j != molecule.atoms.end(); j ++)
            {
                Atom other = (*j);
                
                if (it != j && chem.AtomsCanBond(atom, other))
                {
                    auto bond = make_tuple(it-molecule.atoms.begin(), j-molecule.atoms.begin());
                    molecule.bonds.push_back(bond);
                }
            }
        }
        
        return molecule;
    }
    
    Molecule ParseData(string file)
    {
        ifstream stream(file);
        return ParseData(stream);
    }
    PeriodicTable chem;
};

//View Details for ModelView Matrix
struct View
{
    View(double x, double y, double z, double tx, double ty, double tz)
    : eyeX(x), eyeY(y), eyeZ(z), targetX(tx), targetY(ty), targetZ(tz), xAngle(90), yAngle(270)
    { }
    double eyeX,eyeY,eyeZ,targetX,targetY,targetZ,xAngle,yAngle;
    int width,height;
};

//Generate a rbg value for passed atom
void ColorForType(Symbol type, double* r, double* g, double*b)
    {
        int index = type % 6;
        int row = type / 6;
        if (index == 0 || index == 3 || index == 4) {
            *r = 1.0/(row + 1);
        }
        if (index == 1 || index == 3 || index == 5) {
            *g = 1.0/(row + 1);
        }
        if (index == 2 || index == 4 || index == 5) {
            *b = 1.0/(row + 1);
        }
    }

//render string at passed cordinated with passed rgb value
void RenderString(float x, float y, float z, void *font,string str,float r,float g, float b)
    {  
    char cstr[str.size() + 1];
    strcpy(cstr, str.c_str());
    unsigned char* t = reinterpret_cast<unsigned char *>(cstr);
  glColor3f(r, g, b); 
  glRasterPos3f(x, y, z);
  glutBitmapString(font, t);
    }

//render cylinder by getting the distance and angle between the passed vertices then use it to create a cylinder followed by disc at endpoints
void renderCylinder(float x1, float y1, float z1, float x2,float y2, float z2, float radius,int subdivisions,GLUquadricObj *quadric) 
{
float vx = x2-x1;
float vy = y2-y1;
float vz = z2-z1;

if(vz == 0)
    vz = .0001;

float v = sqrt( vx*vx + vy*vy + vz*vz );
float ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
float rx = -vy*vz;
float ry = vx*vz;
glPushMatrix();

glTranslatef( x1,y1,z1 );
glRotatef(ax, rx, ry, 0.0);
gluQuadricOrientation(quadric,GLU_OUTSIDE);
gluCylinder(quadric, radius, radius, v, subdivisions, 1);

gluQuadricOrientation(quadric,GLU_INSIDE);
gluDisk( quadric, 0.0, radius, subdivisions, 1);
glTranslatef( 0,0,v );

gluQuadricOrientation(quadric,GLU_OUTSIDE);
gluDisk( quadric, 0.0, radius, subdivisions, 1);
glPopMatrix();
}

//render string at passed cordinated with passed rgb value, it change the x value with every new line
void stringPrint(string s,int x,int y,int r,int g, int b)
{
   void* font = GLUT_BITMAP_9_BY_15;
   int qp =0;

     glRasterPos2i(x, y);

    for (string::iterator i = s.begin(); i != s.end(); ++i)
    {
        char c = *i;
        if(c=='\n')
        {
            qp+=20;
            glRasterPos2i(x, y-qp);
        }
        else
        {
            glColor3d(r, g, b);
            glutBitmapCharacter(font, c);
        }
    }
}

//allow printing only if textmode is on
void stringPrinter(string s,int x,int y,int r,int g, int b)
{
    if(textMode){
        stringPrint(s,x,y,r,g,b);
    }
}

//Periodic table
PeriodicTable chem;
//Parser
DataParser parser(chem);
//Initializing View
View view(-22, 0, 0, 0, 0, 0);


//Render the atom and bonds
void drawReaction(){

Molecule molecule[3] = {parser.ParseData(mole[moleI].a),parser.ParseData(mole[moleI].b),parser.ParseData(mole[moleI].c)};
int pk =presentation;
for(int i=0,k=-1; i<3; i++,k++)
        {
            if(moleI > 4 && presentation != -2) //check for presentation mode
            {
                presentation = -2;
                glutPostRedisplay();
            }
            if(pk == 0)
            {
                stringPrinter(projectInfo,400,450,0.0,1.0,1.0);
            }
            if(pk   == -1 || pk > 0)
            {
                for (auto it = molecule[i].atoms.begin(); it != molecule[i].atoms.end(); it ++)
        {
            if(atomMode) //check if atom is to be shown
            {
                Atom& atom = (*it);
                auto it1 = find_if(chem.table.cbegin(), chem.table.cend(), [&atom](AtomDetail x) {
            return x.symbol == atom.symbol;
        }); 

        float radius = (*it1).radius;

            double r = 0.0, g = 0.0, b = 0.0;

            glStencilFunc(GL_ALWAYS,atom.symbol, -1);
            // ColorForType(atom.symbol, &r, &g, &b);
            glColor3d(1, 1, 1);
            
            glPushMatrix();
            glEnable(GL_TEXTURE_2D);
  	        glBindTexture(GL_TEXTURE_2D, texture[atom.symbol]);
            
            GLUquadric *quad;
            quad = gluNewQuadric();
            gluQuadricTexture(quad, 40);
            glTranslatef(atom.x- (10*k),atom.y,atom.z);
            gluSphere(quad, radius*atomSize, 100, 100);
            
            glDisable(GL_TEXTURE_2D);
            glPopMatrix();
            }
        
        if (bondMode) //check if bond is to be shown
        {
            glColor3d(0.9, 1.0, 0.7);
            for (auto it = molecule[i].bonds.begin(); it != molecule[i].bonds.end(); it ++)
            {
                Atom atom = molecule[i].atoms[get<0>(*it)];
                Atom other = molecule[i].atoms[get<1>(*it)];
                
                glStencilFunc(GL_ALWAYS, 0, -1);
                
                GLUquadricObj *quadric=gluNewQuadric();
                gluQuadricNormals(quadric, GLU_SMOOTH);
                renderCylinder(atom.x- (10*k),atom.y,atom.z,other.x- (10*k),other.y,other.z,0.1,40,quadric);
                
                gluDeleteQuadric(quadric);
            }
        }

        if(i==0)
        {
                glColor3d(0.0, 1.0, 0.0);
                
                GLUquadricObj *quadric=gluNewQuadric();
                gluQuadricNormals(quadric, GLU_SMOOTH);
                renderCylinder(4,0,0,5,0,0,0.1,40,quadric);
                gluDeleteQuadric(quadric);
                
                quadric=gluNewQuadric();
                gluQuadricNormals(quadric, GLU_SMOOTH);
                renderCylinder(4.5,-0.5,0,4.5,0.5,0,0.1,40,quadric);
                gluDeleteQuadric(quadric);
        }
        else
        if(i==1)
        {
                glColor3d(1.0, 0.0, 0.0);
                
                GLUquadricObj *quadric=gluNewQuadric();
                gluQuadricNormals(quadric, GLU_SMOOTH);
                renderCylinder(-4.5,0.25,0,-5.5,0.25,0,0.1,40,quadric);
                gluDeleteQuadric(quadric);
                
                quadric=gluNewQuadric();
                gluQuadricNormals(quadric, GLU_SMOOTH);
                renderCylinder(-4.5,-0.25,0,-5.5,-0.25,0,0.1,40,quadric);
                gluDeleteQuadric(quadric);
        }
           
    }
    if(pk != -1)
    {   
    pk --;
    }
    if(pk==1 && k==1)
    {
        presentation=1;
        moleI++;
        glutPostRedisplay();
    }
            }
    }
}

//rendering flooe
void drawFloor(){
if(floorMode)
{

    glBegin(GL_QUADS);
    glVertex3fv(fVert[0]);
    glVertex3fv(fVert[1]);
    glVertex3fv(fVert[2]);
    glVertex3fv(fVert[3]);
    glEnd();
}
}

//main display function
void display(){
    //checking for menu clicks
    if(value == 20){
      atomMode = false;
        value = 0;
        createMenu(); 
    }
    if(value == 21){
      atomMode = true;
        value =0;
        createMenu();
    }
     if(value == 30){
       bondMode = false;
        value = 0;
        createMenu(); 
    }
    if(value == 31){
       bondMode = true;
        value =0;
        createMenu();
    }
    if(value == 40){
        atomSize += 0.1;
        value = 0;
        createMenu(); 
    }
    if(value == 41){
        atomSize -= 0.1;
        value =0;
        createMenu();
    }
     if(value == 50){
        presentation = -1;
        moleI=0;
        value = 0;
        createMenu(); 
    }
    if(value == 51){
        presentation = 0;
        moleI =0;
        value =0;
        createMenu();
    }
    if(value == 60){
        textMode = 0;
        value = 0;
        createMenu(); 
    }
    if(value == 61){
        textMode = 1;
        value = 0;
        createMenu(); 
    }
    if(value == 70){
        floorMode = 0;
        value = 0;
        createMenu(); 
    }
    if(value == 71){
        floorMode = 1;
        value = 0;
        createMenu(); 
    }
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT| GL_STENCIL_BUFFER_BIT);
    glLoadIdentity();
    
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

    gluLookAt(view.eyeX, view.eyeY, view.eyeZ, view.targetX, view.targetY, view.targetZ, 0, 0,1);
    glRotatef(view.xAngle, 0.0f, 0.0f, 1.0f);
    glRotatef(view.yAngle, 1.0f, 0.0f, 0.0f);
    glPushMatrix();
    
    glColor4f(0.3, 0.3, 0.3, 1.0);
    drawFloor();
    glPopMatrix();

    glPushMatrix();
    drawReaction();
    glPopMatrix();
   
    
    glMatrixMode(GL_PROJECTION);
    
    glPushMatrix();
    glLoadIdentity();

    gluOrtho2D(0.0, 1080, 0.0, 700);
    glMatrixMode(GL_MODELVIEW);
    
    glPushMatrix();
    glLoadIdentity();
    if(presentation != -2 && (moleI > 0 || presentation !=0))
    {
        if(presentation == 0)
        {
            stringPrinter(mole[moleI-1].name,850,650,1.0,0.0,1.0);
        }
        else
        {
            stringPrinter(mole[moleI].name,850,650,1.0,0.0,1.0);
        }
    stringPrinter(atomInfoS,10,130,0.0,1.0,0.0);
    stringPrinter(projectInfo,10,660,0.0,1.0,1.0);
    stringPrinter(Instruction,700,120,1.0,1.0,0.0);
    }
    if(presentation == 0 && moleI == 0)
    {
        stringPrinter(projectInfo,430,450,0.0,1.0,1.0);
    }
    if(presentation == -2)
    {
        stringPrinter("-->\nSlide Finished click to exit presentation mode",430,450,0.0,1.0,1.0);
    }

    glMatrixMode(GL_PROJECTION); 

    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);

    glPopMatrix();

    glutSwapBuffers();

}

//Called on mouse click and get the stencil index of object on which it is clicked
void getObj(int button, int state, int x, int y){
    if(state != GLUT_DOWN) return;
    if(start ==0)
    {
            presentation = -1;
            start = -1;
            createMenu();
    }else
    {
         int w_height = glutGet(GLUT_WINDOW_HEIGHT);

        GLuint index;

        glReadPixels(x, w_height - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &index);
        atomInfoS = "Atomic Number : "+ to_string(index) +"\n"+AtomInfo[index];
    if(presentation == -2)
    {
        presentation = -1;
        moleI=0;
        createMenu();
    }

         if(presentation != -1 && index == 0 && presentation != -2)
    {
        presentation += 1;
    }
    }
       
    
        glutPostRedisplay(); 
}

//called on window size change
void reshape(int w, int h){
    if(h == 0)    h = 1;
    
    float ratio = w* 1.0 / h;
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0, ratio, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

//special keyboard
void specialKeyboard(int key, int xx, int yy){
    switch(key){
        case(GLUT_KEY_F1):
            break;
            
        case(GLUT_KEY_F2):
            break;
            
        case(GLUT_KEY_UP):
            break;
            
        default: break;
    }
    glutPostRedisplay();  
}

int lastx=0,lasty=0;

//normal key handler
void NormalKeyHandler (unsigned char key, int x, int y)
{
    if (key == 'z')
    {
        if(lastx==0)
        {
            lastx=x;lasty=y;
        }
        else
        {
            view.xAngle += (lastx -  x)*0.4;
            view.yAngle += (lasty - y)*0.4;
            lastx = x; lasty=y;
        }
    }
    if(presentation == -1 )
    {
if (key == 'w')
    {
            view.eyeX += 1;
    }
    if (key == 's')
    {
            view.eyeX -= 1;
    }
    if (key == 'a')
    {
            view.eyeY += 1;
    }
    if (key == 'd')
    {
            view.eyeY -= 1;
    }
    if (key == 'q')
    {
            view.eyeZ += 1;
    }
    if (key == 'e')
    {
            view.eyeZ -= 1;
    }
    if(key=='x')
    {
        lastx=0,lasty=0;
    }
    if(key=='n')
    {
        moleI = (moleI+1)%5;
    }
    if(key=='p')
    {
        moleI = (moleI+4)%5;
    }
    }
    
    glutPostRedisplay(); 
}

//main function
int main(int argc, char** argv)
{
    getRxns();
    detail();

    glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_STENCIL | GLUT_MULTISAMPLE);
    glutInitWindowSize (1080, 700);
    glutInitWindowPosition (50, 50);
    glutCreateWindow ("Computer Graphics Project");
    //texture array
    texture[1] = loadTextureFromFile("./bmp/hydrogen.bmp");
    texture[6] = loadTextureFromFile("./bmp/carbon.bmp");
    texture[7] = loadTextureFromFile("./bmp/nitrogen.bmp");
    texture[8] = loadTextureFromFile("./bmp/oxygen.bmp");
    texture[16] = loadTextureFromFile("./bmp/sulfur.bmp");
    texture[17] = loadTextureFromFile("./bmp/chlorine.bmp");
    texture[35] = loadTextureFromFile("./bmp/bromine.bmp");
    
    createMenu();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(getObj);
    glutKeyboardFunc (NormalKeyHandler);
    glutSpecialFunc(specialKeyboard);

    glutMainLoop();

    return 0;
}
