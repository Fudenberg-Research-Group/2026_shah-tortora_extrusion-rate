import numpy as np


def compute_LEF_pos(extrusion_engine, n_tot,
                    trajectory_length, dummy_steps,
                    LEF_lifetime, LEF_separation, LEF_stall, TAD_size,
                    **kwargs):
    """LEF dynamics computation"""

    LEF_num = n_tot // LEF_separation
    
    birth_array = np.zeros(n_tot, dtype=np.double) + 0.1
    pause_array = np.zeros(n_tot, dtype=np.double)
        
    death_array = np.zeros(n_tot, dtype=np.double) + 1./LEF_lifetime
    stall_death_array = np.zeros(n_tot, dtype=np.double) + 1./LEF_lifetime
    
    stall_list = np.arange(0, n_tot, TAD_size)
    stall_left_array = np.zeros(n_tot, dtype=np.double)
    stall_right_array = np.zeros(n_tot, dtype=np.double)
    
    for i in stall_list:
        stall_left_array[i] = LEF_stall
        stall_right_array[i] = LEF_stall
        
    translocator = extrusion_engine(birth_array, death_array,
                                    stall_left_array, stall_right_array,
                                    pause_array, stall_death_array,
                                    LEF_num, **kwargs)
    
    translocator.steps(dummy_steps)

    LEF_positions = np.zeros((trajectory_length, LEF_num, 2), dtype=int)

    for step in range(trajectory_length):
        translocator.steps(1)
            
        LEF_positions[step] = np.asarray(translocator.getLEFs()).T

    return LEF_positions
    
    
class LEFTranslocatorDirectional():
    
    def __init__(self, emissionProb, deathProb, stallProbLeft, stallProbRight, pauseProb, stallFalloffProb,  numLEF):
        emissionProb[0] = 0
        emissionProb[len(emissionProb)-1] = 0
        
        emissionProb[stallProbLeft > 0.9] = 0        
        emissionProb[stallProbRight > 0.9] = 0        
        
        self.N = len(emissionProb)
        self.M = numLEF
        self.emission = emissionProb
        self.stallLeft = stallProbLeft
        self.stallRight = stallProbRight
        self.falloff = deathProb
        self.pause = pauseProb
        cumem = np.cumsum(emissionProb, dtype=np.double)
        cumem /= cumem[-1]
        self.cumEmission = cumem
        self.LEFs1 = np.zeros((self.M), int)
        self.LEFs2 = np.zeros((self.M), int)
        self.stalled1 = np.zeros(self.M, int)
        self.stalled2 = np.zeros(self.M, int)
        self.occupied = np.zeros(self.N, int)
        self.stallFalloff = stallFalloffProb
        self.occupied[0] = 1
        self.occupied[self.N - 1] = 1
        self.maxss = 1000000
        self.curss = 99999999

        for ind in range(self.M):
            self.birth(ind)


    def birth(self, ind):
    
        while True:
        
            pos = self.getss()
            
            if pos >= self.N - 1:
                print("bad value", pos, self.cumEmission[len(self.cumEmission)-1])
                continue
                
            if pos <= 0:
                print("bad value", pos, self.cumEmission[0])
                continue 
 
            if self.occupied[pos] == 1:
                continue
            
            self.LEFs1[ind] = pos
            self.LEFs2[ind] = pos
            
            self.occupied[pos] = 1
            
            if (pos < (self.N - 3)) and (self.occupied[pos+1] == 0):
                if np.random.random() > 0.5:
                    self.LEFs2[ind] = pos + 1
                    self.occupied[pos+1] = 1
            
            return


    def death(self):
    
        for i in range(self.M):
        
            if self.stalled1[i] == 0:
                falloff1 = self.falloff[self.LEFs1[i]]
            else: 
                falloff1 = self.stallFalloff[self.LEFs1[i]]
                
            if self.stalled2[i] == 0:
                falloff2 = self.falloff[self.LEFs2[i]]
            else:
                falloff2 = self.stallFalloff[self.LEFs2[i]]              
            
            falloff = max(falloff1, falloff2)
            
            if np.random.random() < falloff:
                self.occupied[self.LEFs1[i]] = 0
                self.occupied[self.LEFs2[i]] = 0
                
                self.stalled1[i] = 0
                self.stalled2[i] = 0
                
                self.birth(i)
    
    
    def getss(self):
    
        if self.curss >= self.maxss - 1:
            foundArray = np.array(np.searchsorted(self.cumEmission, np.random.random(self.maxss)), dtype = np.longlong)
            self.ssarray = foundArray

            self.curss = -1
        
        self.curss += 1
        
        return self.ssarray[self.curss]
        

    def step(self):
        for i in range(self.M):            
            stall1 = self.stallLeft[self.LEFs1[i]]
            stall2 = self.stallRight[self.LEFs2[i]]
                                    
            if np.random.random() < stall1:
                self.stalled1[i] = 1
                
            if np.random.random() < stall2:
                self.stalled2[i] = 1
                         
            cur1 = self.LEFs1[i]
            cur2 = self.LEFs2[i]
            
            if self.stalled1[i] == 0: 
                if self.occupied[cur1-1] == 0:
                    pause1 = self.pause[self.LEFs1[i]]
                    
                    if np.random.random() > pause1:
                        self.occupied[cur1 - 1] = 1
                        self.occupied[cur1] = 0
                        
                        self.LEFs1[i] = cur1 - 1
                        
            if self.stalled2[i] == 0:                
                if self.occupied[cur2 + 1] == 0:                    
                    pause2 = self.pause[self.LEFs2[i]]
                    
                    if np.random.random() > pause2:
                        self.occupied[cur2 + 1] = 1
                        self.occupied[cur2] = 0
                        
                        self.LEFs2[i] = cur2 + 1
        
        
    def steps(self,N):
        for i in range(N):
            self.death()
            self.step()
            
    def getOccupied(self):
        return np.array(self.occupied)
    
    def getLEFs(self):
        return np.array(self.LEFs1), np.array(self.LEFs2)
        
    def updateMap(self, cmap):
        cmap[self.LEFs1, self.LEFs2] += 1
        cmap[self.LEFs2, self.LEFs1] += 1

    def updatePos(self, pos, ind):
        pos[ind, self.LEFs1] = 1
        pos[ind, self.LEFs2] = 1
