



Ez=sqrt(OmegaM*(1+z_)^3+OmegaR*(1+z_)^4+OmegaK*(1+z_)^2+OmegaL*(1+z_)^(3*(1+w)))
SimpsonEz<-z/6*(1/Ez(OmegaM=1,OmegaL=0,w=0,z_=z)+4/Ez(OmegaM=1,OmegaL=0,w=0,z_=(z/2))+1/Ez(OmegaM=1,OmegaL=0,w=0,z_=0))

