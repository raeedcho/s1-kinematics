'''
    # To do in neuromechanic engine:
    # Do preprocessing of nm file : case-sensitive tags
    #                             : _ cannot lead or trail tag name
    
    x = Tree('model')  # Returns the model tree of a loaded Neuromechanic File
    
    Get a RigidBody by it's name attribute:
        x.NeuromechanicFile.Bodies.RigidBody['Name=?']
    
    Get a DegreeOfFreedom position by it's name attribute:
        x.NeuromechanicFile.Bodies.RigidBody.DegreeOfFreedom['Name=?'].State.x_[0]
    
    Get a motor neuron for each muscle
        [x.NeuromechanicFile.Neurons.Neuron['Name='+mn] for mn in x.NeuromechanicFile.Muscles.Muscle.MotorNeuron.Name.x_]
    
'''

import neuromechanic as nm
import numpy as np
from sys import float_info
import math
import collections
import itertools

NODE_ATTRIBUTE=1    # For creating an attribute
NODE_ELEMENT=2      # For creating an element
NODE_COMMENT=3      # For creating a comment
machineepsilon = float_info.epsilon # machine epsilon

class nmerror(Exception):
    nm.system_log(str(Exception))
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def flatten(l):
    return sum([flatten(e) if isinstance(e, list) else [e] for e in l],[])

def xmlstr(x):
    if type(x)==type(np.array([])):
        outvar = str(x.tolist()).replace('[','').replace(']','').replace(',',' ')
    else:
        outvar = str(x)
    return outvar        
            
class TreeList(list):
    
    def __getitem__(self, item):
        
        if isinstance(item,str):
            
            # Get nodes by attribute or element value!
            try:
                name,value = item.split('=');
                outvar = TreeList([])
                for x in self:
                    try:
                        if str(x.__getattr__(name).x_).lower()==value.lower():
                           outvar.append(x) 
                    except:
                        stmt = False
                return outvar
            except:
                raise nmerror('Invalid search criteria: '+item)
            
        else:
            outvar = super(TreeList, self).__getitem__(item)
            if type(outvar)==type([]):
                outvar = TreeList(outvar)
            return outvar
    
    def __getattr__(self, item):
        # If what you are looking for is NOT an endpoint then flatten it and 
        # wrap it as a TreeList
        if len(self)>0:
            if isinstance(self[0].__getattr__(item),TreeList):
                return TreeList(flatten([x.__getattr__(item) if x else None for x in self]))
            elif (len(self)==1):
                return self[0].__getattr__(item)
            else:
                return [x.__getattr__(item) for x in self if x]
        else:
            return TreeList([]);  # Return empty treelist
    
    def addElement_(self,name,value=None):
        [x.addElement_(name,value) for x in self if x]
    def addAttribute_(self,name,value):
        [x.addAttribute_(name,value) for x in self if x]
    def addComment_(self,value):
        [x.addComment_(value) for x in self if x]
        
    def __setattr__(self,name,value):
        if name in ['x_']:
            # Update object in neuromechanic memory
            [x.__setattr__(name,value) for x in self if x]
        else:
            # default setattr
            super().__setattr__(name, value)

class Tree():
    
    def __init__(self,x=[]):
        if isinstance(x,str):
            x = nm.tree_gettreebyname(x);
        
        if x==[]:
            pass # empty class
            
        self.__nodep = x;
        self.__nodechildren = {};
        self.__nodeattributes = {};
        self.__nodecomments = {};
        self.__tag = nm.tree_gettag(x).lower();
        for ii in range(nm.tree_getlength(x)):
            p = nm.tree_getchildi(x,ii);
            tag = nm.tree_gettag(p).lower();    # want everything to be case insensitive
            typ = nm.tree_getnodetype(p);
            
            # These are all lists!
            if NODE_ELEMENT == typ:
                if tag.lower() not in self.__nodechildren:
                    self.__nodechildren[tag.lower()] = TreeList([]);
                xTree = Tree(p)
                self.__nodechildren[tag.lower()].append(Tree(p))
            elif NODE_ATTRIBUTE == typ:
                self.__nodeattributes[tag] = Tree(p);
            elif NODE_COMMENT == typ:
                self.__nodecomments[tag] = Tree(p);
    
    def __getitem__(self,item):
        if isinstance(item,str):
            return self.__getattr__(item)
        
    def addElement_(self,name,value=None):
        if name.lower() not in self.__nodechildren:
            self.__nodechildren[name.lower()] = TreeList([]);
        p = nm.tree_createnode(self.__nodep,NODE_ELEMENT,name)
        self.__nodechildren[name.lower()].append(Tree(p))
        if type(value)==type(np.array([])):
            nm.tree_setvalue(p,value.tolist())
        elif value is not None:
            nm.tree_setvalue(p,value)
        
    def addAttribute_(self,name,value):
        p = nm.tree_createnode(self.__nodep,NODE_ATTRIBUTE,name)
        self.__nodeattributes[name.lower()] = Tree(p)
        if type(value)==type(np.array([])):
            nm.tree_setvalue(p,value.tolist())
        elif value is not None:
            err = nm.tree_setvalue(p,value)
        
    def addComment_(self,value):
        p = nm.tree_createnode(self.__nodep,NODE_COMMENT,value)
        self.__nodecomments[value] = Tree(p)
        
    def __setattr__(self,name,value):
        if name=='x_':
            # Update object in neuromechanic memory
            if type(value)==type(np.array([])):
                nm.tree_setvalue(self._Tree__nodep,value.tolist())
            else:
                nm.tree_setvalue(self._Tree__nodep,value)
        else:
            # default setattr
            super().__setattr__(name, value)
            
    def __getattr__(self,itemc):
        
        item = itemc.lower()
        if item=='x_':
            outvar = nm.tree_getvalue(self.__nodep)
            if type(outvar)==type([]):
                outvar = np.array(outvar)
            # if this is a list make it a numpy array
            return outvar
        elif item=='p_':
            return self.__nodep
        elif item=='c_':
            return self.__nodechildren;
        elif item=='a_':
            return self.__nodeattributes;
            
        elif item=='tag_':
            return self.__tag
        elif item=='xs_':
            outvar = self.x_
            return xmlstr(outvar)
        elif item=='print_':
            strng = '<'+self.__tag;
            for x in self.__nodeattributes:
                # The node attributes dictionary is keyed with the attribute name and the value of the dictionary is a TreeList of length 1 element which is the Tree of the attribute itself;
                strng += ' '+x+'='+self.__nodeattributes[x].xs_
            strng += '>';
            for x in self.__nodechildren:
                z = self.__nodechildren[x].print_;
                for zz in z:
                    strng += zz
            value = self.x_
            if value is not None:
                strng += self.xs_
            
            strng += '</'+self.__tag+'>'
            
            return strng
        elif item=='expointer_':
            return nm.tree_getnodeexternalpointer(self.__nodep)
        elif item=='enum1_':
            return nm.tree_getnodeenum1(self.__nodep)
        elif item=='list_':
            nm.system_log('Attributes:')
            for x in self.__nodeattributes:
                # The node attributes dictionary is keyed with the attribute name and the value of the dictionary is a TreeList of length 1 element which is the Tree of the attribute itself;
                nm.system_log('  '+x+'='+self.__nodeattributes[x].xs_)
            nm.system_log('Children:')
            for x in self.__nodechildren:
                nm.system_log('  '+x+' ('+str(len(self.__nodechildren[x]))+')')
            nm.system_log('Value:')
            nm.system_log('  '+self.xs_)
            nm.system_log('  ')
            
            
            # if this is a list make it a numpy array
            return None
        elif item in self.__nodechildren:
            return self.__nodechildren[item];
        elif item in self.__nodeattributes:
            return self.__nodeattributes[item];
            # return nm.tree_getvalue(self.__nodeattributes[item].__nodep);
        elif item in self.__nodecomments:
            return nm.tree_getvalue(self.__nodecomments[item].__nodep);
        else:
            return TreeList([]);  # Return empty treelist
                
class Build():

    def __init__(self,x=[]):
        # Create a new 
        self.__x = Tree('new')
        self.__x.addElement_('NeuromechanicFile')

    def addRigidBody(self,Name,Parent='ground',Inertia=[.1,.1,.1,0,0,0],\
            FrameLocation=[0,0,0],CenterOfMass=[0,0,0],Mass=None,\
            InertialEllipsoid=[1,1,1],R=np.eye(3),Density=None):
        ''' Add a RigidBody to the structure (x). Inputs are the name of the body (Name) and its 
        parent (Parent), the location (FrameLocation) of the bodies coordinate system in 
        its parents coordinate system, the location of the center of mass (CenterOfMass)
        The dimensions of an equivalent ellipsoid ellipsoids principal axes (Idims). And
        a transformation matrix (R)
        '''
        
        # Convert from principal to semi-principal axis lengths
        if (Mass is None) and (Density is not None):
            Idims = np.array(InertialEllipsoid)
            Idims = Idims/2 
            
            # mass = density*(volume of an ellipsoid)
            Mass = 4*np.pi*density*Idims[0]*Idims[1]*Idims[2]/3 
            
            # Get the inertia from the dimensions of the equivalent inertial ellipsoid
            I = Mass/5*np.array([[Idims[1]*Idims[1]+Idims[2]*Idims[2],0,0],
                            [0,Idims[0]*Idims[0]+Idims[2]*Idims[2],0],
                            [0,0,Idims[1]*Idims[1]+Idims[0]*Idims[0]]]);
            I = np.dot(R,np.dot(I,R.T))
        
            # This is the form that neuromechanic wants to see inertia in
            Inertia = np.array([I[0][0],I[1][1],I[2][2],I[0][1],I[1][2],I[0][2]])
        elif Mass is None:
            raise nmerror('Either "Mass" or "Density" must be specified')
        
        
        # initialize Bodies if necessary
        if not 'Bodies'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Bodies')
            
        # Add the body
        self.__x.NeuromechanicFile.Bodies.addElement_('RigidBody')
        B = self.__x.NeuromechanicFile.Bodies.RigidBody[-1]
        B.addAttribute_('Name',Name)
        B.addElement_('FrameLocation',FrameLocation)
        B.FrameLocation.addAttribute_('Parent',Parent)
        B.addElement_('Mass',Mass)
        B.addElement_('Inertia',Inertia)
        B.addElement_('CenterOfMass',CenterOfMass)
        return self
        
    def addDoF(self,ind,Name,Type='Rotation',Axis=[1,0,0],Locked='false',State=[0,0]):
        ''' 
        Add a DegreeOfFreedom to the "ind"th RigidBody, Specification includes
        the Name "Name" of the Degree of freedom, the Type (Type="Translation" 
        or "Rotation", the "Axis", whether the degree of freedom is "Locked" and
        the "State" [Position,Velocity]
        '''
        
        # initialize DegreeOfFreedom if necessary
        self.__x.NeuromechanicFile.Bodies.RigidBody[ind].addElement_('DegreeOfFreedom')

        D = self.__x.NeuromechanicFile.Bodies.RigidBody[ind].DegreeOfFreedom[-1]
        D.addAttribute_('Name',Name)
        D.addAttribute_('Type',Type)
        D.addElement_('Locked',Locked)
        D.addElement_('Axis',Axis)
        D.addElement_('State',State)
        return self
        
    def addContact(self,Name,Type,Parent=None,Ellipsoid=None):
        '''
        Add a contact surface (named "Name" of type "Type" on body "Parent") to the structure.
        Parent is not necessary if Type="Polygon". If Type="Ellipsoid" Then a dictionary must
        be passed in. e.g. Ellipsoid={'PrincipalAxes':[.1,.12,.4],'Rotation':[0,0,0],'Origin':[0,0,0]}
        '''
        if not 'Environment'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Environment')
        if not 'Contacts'.lower() in self.__x.NeuromechanicFile.Environment.c_:
            self.__x.NeuromechanicFile.Environment.addElement_('Contacts')
            
        self.__x.NeuromechanicFile.Environment.Contacts.addElement_('Contact')
        c = self.__x.NeuromechanicFile.Environment.Contacts.Contact[-1]
        c.addAttribute_('Name',Name)
        c.addAttribute_('Type',Type)
        if Parent is not None:
            c.addAttribute_('Parent',Parent)
        if Ellipsoid is not None:
            c.addElement_('Ellipsoid')
            for x in Ellipsoid:
                c.Ellipsoid.addElement_(x,Ellipsoid[x])
        
    def addGravity(self,Gravity=[0.,-9.80665,0.]):
        '''Add gravity to the environment'''
        if not 'Environment'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Environment')
        self.__x.NeuromechanicFile.Environment.addElement_('Gravity',Gravity)
        
    def addPerturbation(x,Name,Parent,Type,Wrench=[0.,0.,0.],ApplicationDirection=[1.,0.,0.],
        ApplicationPoint={'Parent':'ground','x':[0.,0.,0.]},ShapeParameters={'Type':'SquarePulse','x':[.1,.1,.2]}):
        ''' 
        Add a perturbation (name=Name) which force is appliced in the "ApplicationDirection"
        direction in the "Parent" frame of reference. The force is applied to "ApplicationPoint['x']"
        in the "ApplicationPoint['Parent']" frame of reference. The "Type" should be either Force or Moment.
        ShapeParameters is a dictionary with type = "SquarePulse|GaussianPulse|SineWave|splineName" and
        x equal to an nx3 matrix for all types except splineName
        '''
        
        if not 'Environment'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Environment')
        if not 'Perturbations'.lower() in self.__x.NeuromechanicFile.Environment.c_:
            self.__x.NeuromechanicFile.Environment.addElement_('Perturbations')
        self.__x.NeuromechanicFile.Environment.Perturbations.addElement_('Perturbation')
        p = self.__x.NeuromechanicFile.Environment.Perturbations.Perturbation[-1]
        
        p.addAttribute_('Name',Name)
        p.addAttribute_('Parent',Parent)
        p.addAttribute_('Type',Type)
        p.addElement_('EqWrench',Wrench)
        p.addElement_('ApplicationDirection',ApplicationDirection)
        p.addElement_('ApplicationPoint',ApplicationPoint['x'])
        p.ApplicationPoint.addAttribute_('Parent',ApplicationPoint['Parent'])
        p.addElement_('ShapeParameters',ShapeParameters['x'])
        p.ShapeParameters.addAttribute_('Parent',ShapeParameters['Type'])
    
    def addPolygon(x,Name,Parent,Vertices,Facets,Origin=[0.,0.,0.],Rotation=[0.,0.,0.]):
        '''
        Add a polygon with "Name" that is defined in the "Parent" RigidBody coordinate system
        with origin at "Origin" and oriented with euler angles in "Rotation". The
        vertex list is in "Vertices" and facet list is in "Facets", 
        '''
        if not 'Resources'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Resources')
        self.__x.NeuromechanicFile.Resources.addElement_('Polygon')
        p = self.__x.NeuromechanicFile.Resources.Polygon[-1]
        
        p.addAttribute_('Name',Name)
        p.addElement_('Origin',Origin)
        p.Origin.addAttribute_('Parent',Parent)
        p.addElement_('Rotation',Rotation)
        p.addElement_('Vertices',Vertices)
        
        if type(Facets[0])==type([]):
            q = len(Facets[0])==4
        else:
            q = len(Facets)==4
            
        if q:
            p.addElement_('Quads',np.array(Facets).astype('b'))
        else:
            p.addElement_('Tris',np.array(Facets).astype('b'))
    
    def addMuscleSplines(self):
        '''
        Add three splines (Active force length curve (afl), passive force length curve (pfl), 
        and active force velocity curve (afv)) to the Functions section
        '''
        
        afl = [[-5.,0.],[0.,0.],[4.01E-1,0.],[4.02E-1,0.],[4.035E-1,0.],
            [5.2725E-1,2.2667E-1],[6.2875E-1,6.3667E-1],[7.1875E-1,8.5667E-1],
            [8.6125E-1,9.5E-1],[1.045,9.9333E-1],[1.2175,7.7E-1],[1.4387,2.4667E-1],
            [1.6187,0.],[1.62,0.],[1.621,0.],[2.2,0.],[5.,0.]];

        pfl = [[-5.,0.],[9.98E-1,0.],[9.99E-1,0.],[1.,0.],[1.1,3.5E-2],
            [1.2,1.2E-1],[1.3,2.6E-1],[1.4,5.5E-1],[1.5,1.17],[1.6,2.],[1.601,2.],[1.602,2.],
            [5.,2.]]
            
        afv = [[-1., 0.],[-9.5E-1,1.0417E-2],[-9.E-1,2.1739E-2],
            [-8.5E-1,3.4091E-2],[-8.E-1,4.7619E-2],[-7.5E-1,6.25E-2],[-7.E-1,7.8947E-2],
            [-6.5E-1,9.7222E-2],[-6.E-1,1.1765E-1],[-5.5E-1,1.4062E-1],[-5.E-1,1.6667E-1],
            [-4.5E-1,1.9643E-1],[-4.E-1,2.3077E-1],[-3.5E-1,2.7083E-1],[-3.E-1,3.1818E-1],
            [-2.5E-1,3.75E-1],[-2.E-1,4.4444E-1],[-1.5E-1,5.3125E-1],[-1.E-1,6.4286E-1],
            [-5.E-2,7.9167E-1],[0.,1.],[5.E-2,1.482],[1.E-1,1.6016],[1.5E-1,1.6558],
            [2.E-1,1.6867],[2.5E-1,1.7068],[3.E-1,1.7208],[3.5E-1,1.7311],[4.E-1,1.7391],
            [4.5E-1,1.7454],[5.E-1,1.7505],[5.5E-1,1.7547],[6.E-1,1.7583],[6.5E-1,1.7614],
            [7.E-1,1.764],[7.5E-1,1.7663],[8.E-1,1.7683],[8.5E-1,1.7701],[9.E-1,1.7717],
            [9.5E-1,1.7732],[1.,1.7745]];
            
        if not 'Functions'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Functions')
            
        self.__x.NeuromechanicFile.Functions.addElement_('Spline');
        s = self.__x.NeuromechanicFile.Functions.Spline[-1]
        s.addAttribute_('Name','afl')
        s.addAttribute_('Type','NaturalCubic')
        s.addElement_('ConstraintA',0.)
        s.addElement_('ConstraintB',0.)
        s.addElement_('ControlPoints',afl)
        
        self.__x.NeuromechanicFile.Functions.addElement_('Spline');
        s = self.__x.NeuromechanicFile.Functions.Spline[-1]
        s.addAttribute_('Name','pfl')
        s.addAttribute_('Type','NaturalCubic')
        s.addElement_('ConstraintA',0.)
        s.addElement_('ConstraintB',0.)
        s.addElement_('ControlPoints',pfl)
        
        self.__x.NeuromechanicFile.Functions.addElement_('Spline');
        s = self.__x.NeuromechanicFile.Functions.Spline[-1]
        s.addAttribute_('Name','afv')
        s.addAttribute_('Type','NaturalCubic')
        s.addElement_('ConstraintA',0.)
        s.addElement_('ConstraintB',0.)
        s.addElement_('ControlPoints',afv)
        
    def addNeuron(self,Name,OutputBounds=[0.,1.],EqBounds=[0.01,.95],EqCCost=.0,
        Synapse=[{'Name':'Baseline','Type':'Constant','Value':[0,0]}]):
        ''' Add a Neuron to the structure (x) Inputs are the name of the Neuron (Name), 
        equilibration bounds (EqBounds), linear term of the Equilibration cost function (EqCCost), 
        the output bounds of the motor neuron (OutputBounds), and finally a list of dictionaries 
        containing Synapse inputs [{'Name':...,'Type':...,'Value':...)} ...]
        '''
        # Initialize Neurons and Muscles if necessary
        if not 'Neurons'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Neurons')
            
        self.__x.NeuromechanicFile.Neurons.addElement_('Neuron');
        n = self.__x.NeuromechanicFile.Neurons.Neuron[-1]
            
        n.addAttribute_('Name',Name)
        n.addAttribute_('Type','NeuralNet')
        n.addElement_('OutputBounds',OutputBounds)
        n.addElement_('EqBounds',EqBounds)
        n.addElement_('EqCCost',EqCCost)
        for y in Synapse:
            n.addElement_('Synapse',y['Value'])
            n.Synapse[-1].addAttribute_('Name',y['Name'])
            n.Synapse[-1].addAttribute_('Type',y['Type'])
        
    def addZajacMuscle(self,Name,MusclePath,MotorNeuron,MaxForce=100,OptimalFiberLength=.1,
        TendonSlackLength=.1,Timescale=0.1,PassiveDamping=0.001,PennationAngle=0,
        ActiveForceLength='afl',PassiveForceLength='pfl',ActiveForceVelocity='afv'):
        ''' Add a Zajac style Muscle to the structure (x) and an associated motor neuron. 
        Inputs are the name of the muscle (Name) a list containing dictionaries of the 
        attachment points (MusclePath=[{'Parent':'nameofbody','x':[x,y,z]} ...]). The maximum
        isometric force of the muscle (MaxForce) the optimal fiber length (OptimalFiberLength), 
        the tendon slack length (TendonSlackLength), the timescale (Timescale, for computing 
        vmax), the passive damping (PassiveDamping) the pennation angle (PennationAngle) in 
        radians, the name of the active force length curve (ActiveForceLength), the name of 
        the passive force length curve (PassiveForceLength), and the name of the active force 
        velocity curve (ActiveForceVelocity)
        '''
        
        # Initialize Muscles if necessary
        if not 'Muscles'.lower() in self.__x.NeuromechanicFile.c_:
            self.__x.NeuromechanicFile.addElement_('Muscles')
            
        # Initialize Muscle if necessary
        self.__x.NeuromechanicFile.Muscles.addElement_('Muscle');
        m = self.__x.NeuromechanicFile.Muscles.Muscle[-1]
        
        
        # Build the structure
        m.addAttribute_('Name',Name)
        m.addElement_('MaxForce',MaxForce)
        m.addElement_('MotorNeuron',MotorNeuron)
        m.addElement_('Model')
        m.Model.addAttribute_('Type','Zajac')
        m.Model.addElement_('Timescale',Timescale)
        m.Model.addElement_('FunctionReferences')
        m.Model.FunctionReferences.addElement_('ActiveForceLength',ActiveForceLength)
        m.Model.FunctionReferences.addElement_('PassiveForceLength',PassiveForceLength)
        m.Model.FunctionReferences.addElement_('ActiveForceVelocity',ActiveForceVelocity)
        m.Model.addElement_('OptimalFiberLength',OptimalFiberLength)
        m.Model.addElement_('TendonSlackLength',TendonSlackLength)
        m.Model.addElement_('PennationAngle',PennationAngle)
        m.Model.addElement_('PassiveDamping',PassiveDamping)
        m.Model.addElement_('MusclePath')
        for p in MusclePath:
            m.Model.MusclePath.addElement_('Point',p['x'])
            [m.Model.MusclePath.Point[-1].addAttribute_(x,p[x]) for x in p if x.lower()!='x']
    
    def writeLoad(self,fname):
        # First write to disk
        err = nm.tree_writenodechildren(self.__x.p_,fname)
        nm.tree_destroynode(self.__x.p_)
        if err==0:
            err = nm.file_load(fname)
        return err

def tsladjust(nmtree,tslmin,lftargmax):
    '''
    Adjust all muscles to have a maximum (dimensionless) operating length of lftargmax, 
    with a minimum tendon slack length of tslmin (in units of length e.g. meters)
    '''
    for z in nmtree.NeuromechanicFile.Muscles.Muscle:
        MTL = z.MusculoTendonLength.x_
        LF0 = z.Model.OptimalFiberLength.x_
        TSL = z.Model.TendonSlackLength.x_
        PEN = z.Model.PennationAngle.x_
        mxMTL = MTL;  #max(MTL[:,ii]) # Get the maximum muscle length
        mnMTL = MTL;  #min(MTL[:,ii]) # Get the minimum muscle length
        lfmax = math.sqrt((mxMTL-TSL)**2 + (LF0*math.sin(PEN))**2)
        lfmin = math.sqrt((mnMTL-TSL)**2 + (LF0*math.sin(PEN))**2)
        
        # Solve for TSL to give lfmax/lf0 = 1.05; Minimum TSL = 5mm
        TSL = max(tslmin,mxMTL - math.sqrt(lftargmax**2 - math.sin(PEN)**2)*LF0) 
        z.Model.TendonSlackLength.x_ = TSL
        
        # See what lfmax and lfmin will be
        lfmax = math.sqrt((mxMTL-TSL)**2 + (LF0*math.sin(PEN))**2)
        lfmin = math.sqrt((mnMTL-TSL)**2 + (LF0*math.sin(PEN))**2)
        nm.system_log(z.Name.x_ + ': [LFmax/LF0, LFmin/LF0] = [' + str(lfmax/LF0)+' , '+str(lfmin/LF0) + ']')        
        
def boundingbox(nmtree):
    xmax = [-1.e30,-1.e30,-1.e30]
    xmin = [1.e30,1.e30,1.e30]
    for z in nmtree.NeuromechanicFile.Resources.Polygon:
        b = z.Origin.Parent.enum1_
        if b > 0:  # Make sure that the polygon isn't on the ground body
            fct,V = getMesh(z)
            if not V is None:
                for jj in range(V.shape[0]):
                    x = nm.model_calclocaltoglobal(b,V[jj,0:3].tolist())
                    xmax = [ max(x[kk],xmax[kk]) for kk in range(3)]
                    xmin = [ min(x[kk],xmin[kk]) for kk in range(3)]
    return(xmax,xmin)
    
def print(*arg):
    ''' Prints whatever arg is to the neuromechanic logfile and log window '''
    if (type(arg)==type(' ')):
        nm.system_log(arg)
    else:
        try:
            x = ' '.join(map(str,arg))
        except:
            x = str(arg)
        nm.system_log(x)
        
def crossmat(x):
   return np.array([[0., -x[2], x[1]],\
                    [x[2], 0., -x[0]],\
                    [-x[1], x[0], 0.]])

def pglobal(p):
   return np.array(nm.model_calclocaltoglobal(p.Parent.enum1_,p.x_.tolist()))

def aglobal(nmtree,bod,axis):
   return np.dot(nmtree.NeuromechanicFile.Bodies.RigidBody['name='+bod].RotationMatrix.x_,np.array(axis))
   
def vglobal(p):
   return np.array(nm.model_calcglobalvelocity(p.Parent.enum1_,p.x_.tolist()))
   
def contacts(x):
   # A is a matrix used to calculate the impact of contact forces on center of mass dynamics
   # A*Fctc = [sum(ContactForces) sum(r x ContactForces)]'
   com = np.array(x.NeuromechanicFile.Dynamic.ModelCOM.x_)
   A = []
   
   # Since nodes in active contacts can be added or destroyed during simulation you should recompute
   # the model tree every time it is called
   acp = Tree(x.NeuromechanicFile.Dynamic.ActiveContacts[0]._Tree__nodep)
   for aci in acp.ContactPoint:
      if 1 == aci.IsActive.x_:          # Go through each active contact
         pg = pglobal(aci.Point)        # Get the global coordinates of the contact point
         dcom = [com[ii]-pg[ii] for ii in range(len(pg))] # Find the vector from the contact point to the center of mass of the model
         Ar = np.vstack((np.eye(3),crossmat(dcom))) # Build the matrix A including the cross product of dcom
         if A == []:
            A = Ar
         else:
            A = np.hstack((A,Ar))
   return A
   
class modeldynamics():
   '''
   This class is meant to retrieve dynamic variables from engine as you need them
   M = inertia matrix
   Mi = inverse inertia matrix
   X = State vector
   G = Generalized gravitational torques
   C = Generalized Coriolis... torques
   L = momentum matrix. i.e. L*qdot = [m*vcom ; I*wcom]
   R = Momentarm matrix = transpose(dMuscleLength/dq)
   mass = mass of model
   unl = vector of indices of unlocked dofs
   qdd = Acceleration vector of generalized coordinates
   qd = Velocity vector of generalized coordinates
   Fm = Vector of muscle forces
   Jc = concatenated Jacobians of all kinematic constraints
   JcdotQdot = Time derivative of concatenated constraint endpoint jacobians times velocity of kinematic states
   Jcxdd = the acceleration of constraints as in: Jc*Qdd + JcdotQdot = Jcxdd
   Fc = Vector of concatenated constraint forces of all kinematic constraints
   Jf = concatenated endpoint Jacobians of all contact points
   JfdotQdot = time derivative of concatenated endpoint Jacobians of all contact points times velocity of kinematic states
   Ff = Vector of concatenated ground reaction forces of all contact points
   '''
   def __init__(self,x):
      self._x = x;
      self._dyn = x.NeuromechanicFile.Dynamic;
      self._int = x.NeuromechanicFile.Internal;
      self.nbod = self._int.numBodies.x_       # Retrieve the number of muscles
      self.nmus = self._int.numMuscles.x_       # Retrieve the number of muscles
      self.ndof = self._int.numKinematicDOF.x_  # Retrieve the number of kinematic degrees of freedom
      self._alias = {
        'M':'GenInertia',
        'Mi':'GenInertiaInverse',
        'CoM':'ModelCOM',
        'X':'ModelState',
        'Q':'...',
        'Qd':'...',
        'G':'GravityGenForce',
        'C':'CoriolisGenForce',
        'L':'GenMomentum',
        'R':'MuscleMomentArm',
        'mass':'ModelMass',
        'unl':'UnlockedDOFIndex',
        'Jc':'ConstraintJacobian',
        'JcdotQdot':'ConstraintDJacobianDState',
        'Jcxdd':'ConstraintAcceleration',
        'Jf':'ContactJacobian',
        'Ff':'ContactForce',
        'JfdotQdot':'ContactDJacobianDState',
        'Fm':'MuscleForce',
      }
      
      # This will store dynamic variables that have been retrieved already:
      self._retrieved = {};  
   
   def __getattr__(self,name):
      if name in self._retrieved:  
         return self._retrieved[name]
      elif name in self._alias:
         
         if name=='Mi':
            nm.model_calcinertiainverse()
            z = np.array(getattr(self._dyn,self._alias[name]).x_)
         elif name=='Q':
            z = self.X[0:self.ndof]
         elif name=='Qd':
            z = self.X[self.ndof:2*self.ndof]
         elif name=='unl':
            z = np.array(getattr(self._int,self._alias[name]).x_)-1
         else:   
            z = np.array(getattr(self._dyn,self._alias[name]).x_)
            
         self._retrieved[name] = z
         return z
         
def mergeJacobians(x,dyn):

   ncon = 0
   if (not dyn.JcdotQdot==None):
      ncon = len(dyn.JcdotQdot);
   
   Af = contacts(x)
   if Af==[]:
      nctc = 0
   else:
      nctc = int(Af.shape[1]/3)  # Calculate number of contacts
      
   # Default output is zero
   J = [];
   JdotQdot_xdd = [];
   Fe = [];
   
   if nctc > 0:
      # If there are frictional contacts
      J = dyn.Jf
      JdotQdot_xdd = dyn.JfdotQdot
      Fe = dyn.Ff
      
   if (ncon>0):
      # If there are kinematic constraints enhance J, JdotQdot, and Xnow
      JcRHS = dyn.JcdotQdot - dyn.Jcxdd # Retrieve the current contstraint Jacobians
      if nctc == 0:
         J = dyn.Jc
         JdotQdot_xdd = JcRHS
         Fe = dyn.Fc
      else:
         J = np.vstack((J,dyn.Jc))
         JdotQdot_xdd = np.hstack((JdotQdot_xdd,JcRHS))
         Fe = np.hstack((Fe,dyn.Fc))

   return (nctc,ncon,J,JdotQdot_xdd,Fe,Af)
    
def ischildof(nmtree,childname,parentname):
    ch = nmtree.NeuromechanicFile.Bodies.RigidBody['Name='+childname]
    while ch is not None:
        parnt = ch.FrameLocation.Parent.x_
        if parnt.lower()==parentname.lower():
            return True
        elif parnt.lower()=='ground':
            return False
        else:
            ch = nmtree.NeuromechanicFile.Bodies.RigidBody['Name='+parnt]

def delaysignal(t,D,x=[]):
    ''' 
    Generates & uses a circular buffer for delaying signal "x"
    t the current time
    x the value of the signal at time t
    y the interpolated value of the signal at time t-tdelay

    D['delay'] the delay of signal x in time 
    D['interval'] Time is discretized and interpolated with this value, D['delay'] should be
                 divisible by D['interval']
    D['buffer'] is the actual buffer and should be initialized to []. The first
                  time that delaysignal is called the buffer will be populated with x
                  mimicing a steady-state of x before the simulation begins.
                  The buffer looks something like this
                  D['buffer']  = [S_-1, S_0, S_1, S_2, S_3, S_4, S_5, ..., S_n]
    D['timen'] will be populated with the initial value of time the first time that
                delaysignal is called. Is the time at which event S_n happened
                
    Here's some code to test it, just delays the time signal itself
    dly = {}
    dly['delay'] = .1
    dly['interval'] = .001
    dly['buffer'] = []
    dly['timen'] = 0.

    # Build t which simulates the time intervals in rk45
    tsub = [.25,.5,.75,1.001,.875]
    tspan = [0]
    n = int(2.*np.floor(dly['delay']/dly['interval']))
    for ii in range(n):
        for ts in tsub:
            tspan.append((ii+ts)*dly['interval'])
            
    delaysignal(0,dly,[0.])
    for t in tspan:
        td = delaysignal(t,dly) # Get delayed value
        if t>=dly['delay']:
            print([t,t-td[0]]) # As long as t >= dly['delay'] this should print a number very close to the delay
        delaysignal(t,dly,[t])  # Potentially add delay to buffer
    
                
    '''
    
    # Initialize the buffer it is not already initialized
    if D['buffer']==[]:
        np1 = int(np.floor(D['delay']/D['interval']+.5)) + 2
        if np1>2 and len(x)!=0:  #  No delay or nothing to populate the buffer with
            D['buffer'] = collections.deque(maxlen=np1)
            D['timen'] = t
            for ii in range(np1):
                D['buffer'].append(np.array(x))
        return x
    
    # Retrieve delayed value
    rmndr = (t - D['timen'])/D['interval']
    if rmndr < 0. and rmndr > -1.+machineepsilon:
        prcnt = rmndr+1
        y = prcnt*D['buffer'][1] + (1.-prcnt)*D['buffer'][0]
    elif rmndr >= 0.:
        n = int(np.floor(rmndr))
        prcnt = rmndr-n
        y = prcnt*D['buffer'][n+2] + (1.-prcnt)*D['buffer'][n+1]
    else:
        y = D['buffer'][0]
        #raise nmerror('Problem with delay buffer t='+str(t)+' tn='+str(D['timen']))
        
    # Store x in buffer
    if len(x) > 0 and t+machineepsilon > D['timen'] + D['interval']: 
        prcnt = (t - (D['timen'] + D['interval']))/D['interval']
        xint = (1.-prcnt)*np.array(x) + prcnt*D['buffer'][-1]
        D['timen'] = D['timen'] + D['interval']
        D['buffer'].append(xint)
            
    return y
    
def getMesh(ply,fct=None,vrt=None):

    if len(ply.Vertices) == 0:
        if len(ply.File)==0:
            raise nmerror('This polygon is not properly formatted')
            return
        polyfile = Tree(ply.File.exPointer_).Bones
    else:
        polyfile = ply
    
    if (vrt==None and fct==None):
        vrt = polyfile.Vertices.x_
        
        if len(polyfile.Tris)==0:
            fct = polyfile.Quads.x_
        else:
            fct = polyfile.Tris.x_
        return fct,vrt
    else:
        
        if type(fct[0])==type([]):
            nf = len(fct[0])
        elif type(fct[0])==type(0):
            nf = len(fct)
        else:
            return None,None
            
        if nf==4:
            polyfile.Vertices.x_ = vrt
            polyfile.Quads.x_ = fct
            if len(polyfile.Tris) > 0: # This may be dangerous or wrong
                nme.nm.tree_destroynode(polyfile.Tris._Tree__nodechildren)
        elif nf==3:
            polyfile.Vertices.x_ = vrt
            polyfile.Tris.x_ = fct
            if len(polyfile.Quads) > 0: # This may be dangerous or wrong
                nme.nm.tree_destroynode(polyfile.Quads._Tree__nodechildren)
        return None,None
                
def R2axis(M):

    m00 = M[0, 0]
    m01 = M[0, 1]
    m02 = M[0, 2]
    m10 = M[1, 0]
    m11 = M[1, 1]
    m12 = M[1, 2]
    m20 = M[2, 0]
    m21 = M[2, 1]
    m22 = M[2, 2]
    # symmetric matrix K
    K = np.array([[m00-m11-m22, 0.0,         0.0,         0.0],
                  [m01+m10,     m11-m00-m22, 0.0,         0.0],
                  [m02+m20,     m12+m21,     m22-m00-m11, 0.0],
                  [m21-m12,     m02-m20,     m10-m01,     m00+m11+m22]])
    K /= 3.0
    # quaternion is eigenvector of K that corresponds to largest eigenvalue
    w, V = np.linalg.eigh(K)
    q = V[[3, 0, 1, 2], np.argmax(w)]
    
    # quaternion to Axis
    nrm = np.linalg.norm(q[1:4])
    if nrm<=0.:
        a = np.array([0.,0.,0.])
    elif q[0]<0.:
        a = (np.pi-np.arcsin(nrm))*2*q[1:4]/nrm;
    else:
        a = np.arcsin(nrm)*2*q[1:4]/nrm;
        
    return a

def statePosture(nmtree,name):
    '''
    Set state to a posture. if the posture name is not defined for a given state
    the state will remain as the current value. Right now only for kinematic states
    '''
    X = nmtree.NeuromechanicFile.Dynamic.ModelState.x_ 
    doflist = nmtree.NeuromechanicFile.Bodies.RigidBody.DegreeOfFreedom;
    ndof = nmtree.NeuromechanicFile.Internal.numKinematicDOF.x_;
    
    for dof in doflist:
        xp = dof.Posture['Name='+name].x_;
        if not xp==[]:
            n = dof.DOFIndex.x_;  # DOFS are indexed from 1 in the engine, python likes 0
            X[n] = xp[0]
            X[n+ndof] = xp[1]
    nmtree.NeuromechanicFile.Dynamic.ModelState.x_ = X
    
def exportSTL(node,fname):

    fct,vrt = getMesh(node)
    stlfile = open(fname,'w')
    
    stlfile.write('solid test\n');
    if len(fct.shape)==1:
        # only 1 facet
        nc = fct.shape[0]
        nr = 1;
    else:
        nc = fct.shape[1]
        nr = fct.shape[0]
        
    for jj in range(nr):
        nrm = np.cross(vrt[fct[jj,2]-1,0:3]-vrt[fct[jj,1]-1,0:3],vrt[fct[jj,0]-1,0:3]-vrt[fct[jj,1]-1,0:3])
        nrm = nrm / np.linalg.norm(nrm)
        
        stlfile.write('facet normal '+xmlstr(nrm)+'\n');
        stlfile.write('  outer loop \n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,0]-1,0:3])+'\n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,1]-1,0:3])+'\n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,2]-1,0:3])+'\n');
        stlfile.write('  endloop \n');
        stlfile.write('endfacet \n');
        if nc==3: 
            continue
        if fct[jj,nc-2]==fct[jj,nc-1]:
            continue
            
        nrm = np.cross(vrt[fct[jj,0]-1,0:3]-vrt[fct[jj,3]-1,0:3],vrt[fct[jj,2]-1,0:3]-vrt[fct[jj,3]-1,0:3])
        nrm = nrm / np.linalg.norm(nrm)
        
        stlfile.write('facet normal '+xmlstr(nrm)+'\n');
        stlfile.write('  outer loop \n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,2]-1,0:3])+'\n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,3]-1,0:3])+'\n');
        stlfile.write('    vertex '+xmlstr(vrt[fct[jj,0]-1,0:3])+'\n');
        stlfile.write('  endloop \n');
        stlfile.write('endfacet \n');
    
    stlfile.write('endsolid test\n');
    stlfile.close()

