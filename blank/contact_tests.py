from wall import make_contact_constraint
import numpy as np

def probe_contact(dofs,const_n,constraints,cases):
    for case in cases:
        test_disp = np.zeros((32,1))
        test_disp[dofs[0]]=case['disp'][0]
        test_disp[dofs[1]]=case['disp'][1]
        out = np.matmul(constraints,test_disp)
        print(f"testing disp ({case['disp'][0]},{case['disp'][1]}) on DOF ({dofs}) with constraint {const_n} expecting {case['valid']}")
        assert(out[const_n,0]*case["valid"]>=0)
        print("passed!\n\n")

def contact_gen(bar_1,bar_2):

    const = make_contact_constraint(
        bar_2["pos"],
        bar_1["pos"],
        bar_2["dof"],
        bar_1["dof"],
        8
    )[0,:]
    print(const)
    return const

def ab_contact():
    bar_1 = {
        "pos":np.array([0,0,0,4],dtype=float),
        "dof": [0,1,2,3]
    }
    bar_2 = {
        "pos":np.array([1,2,5,2],dtype=float),
        "dof": [4,5,6,7]
    }
    const = contact_gen(bar_1,bar_2)

    valid_motion = np.zeros((8,1))
    valid_motion[0,0]=1
    valid_motion[2,0]=1
    valid_motion[4,0]=1
    valid_motion[6,0]=1

    assert(np.matmul(const,valid_motion)>=0)

    valid_motion = np.zeros((8,1))
    valid_motion[0,0]=0
    valid_motion[2,0]=0
    valid_motion[4,0]=1
    valid_motion[6,0]=1
    assert(np.matmul(const,valid_motion)>=0)

    invalid_motion = np.zeros((8,1))
    invalid_motion[0,0]=1
    invalid_motion[2,0]=1
    invalid_motion[4,0]=0
    invalid_motion[6,0]=0
    assert(np.matmul(const,invalid_motion)<0)

    invalid_motion = np.zeros((8,1))
    invalid_motion[0,0]=0
    invalid_motion[2,0]=0
    invalid_motion[4,0]=-1
    invalid_motion[6,0]=-1
    assert(np.matmul(const,invalid_motion)<0)
    
def ac_contact():
    bar_1 = {
        "pos":np.array([0,0,0,4],dtype=float),
        "dof": [0,1,2,3]
    }
    bar_2 = {
        "pos":np.array([5,2,1,2],dtype=float),
        "dof": [6,7,4,5]
    }
    const = contact_gen(bar_1,bar_2)

    valid_motion = np.zeros((8,1))
    valid_motion[0,0]=1
    valid_motion[2,0]=1
    valid_motion[4,0]=1
    valid_motion[6,0]=1

    assert(np.matmul(const,valid_motion)>=0)

    valid_motion = np.zeros((8,1))
    valid_motion[0,0]=0
    valid_motion[2,0]=0
    valid_motion[4,0]=1
    valid_motion[6,0]=1
    assert(np.matmul(const,valid_motion)>=0)

    invalid_motion = np.zeros((8,1))
    invalid_motion[0,0]=1
    invalid_motion[2,0]=1
    invalid_motion[4,0]=0
    invalid_motion[6,0]=0
    assert(np.matmul(const,invalid_motion)<0)

    invalid_motion = np.zeros((8,1))
    invalid_motion[0,0]=0
    invalid_motion[2,0]=0
    invalid_motion[4,0]=-1
    invalid_motion[6,0]=-1
    assert(np.matmul(const,invalid_motion)<0)

if __name__=="__main__":
    ab_contact()
    ac_contact()