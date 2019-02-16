# what is radius at water surface for a given depth? assuming vertical bucket
radius<- function(d) 1.5*d/35.668+12.85

# what is volume of water in bucket for a given depth? assuming vertical bucket
volume<- function(d) {
  pi/3*((d+328.38)*(radius(d))^2-328.38*12.85^2)
}

# volume of one sample dip, assuming skirt was FLUSH not high
dip.vol.flush<- function(d) {
  44.683*d-192.97
}
# volume of one sample dip, assuming skirt was HIGH not flush
dip.vol.hi<- function(d) {
  42.8608*d-132.4778
}

# function to calculate new depth after removing a specified sample volume, given initial depth
new.depth<-function(d.start,v.sample){
  v.new<-volume(d.start)-v.sample # calculate new volume as (initial vol.) - (sample vol.)
  return(polyroot(c(-3*v.new/pi,520.04,1.66616,0.0017686))) # I will only need the single positive real root
}

# function to take first root of new.depth fxn, and discard the imaginary part (which should be zero!) 
numeric.depth<-function(d){
  suppressWarnings(as.numeric(d[1])) # don't want to hear about discarded imaginary parts
}

# what is the total volume of a pooled sample, starting from depth=d and taking n dips?
sample.vol.flush<-function(depth,ndips=3){
  d<-depth
  # print(d)
  collected=0
  for(i in 1:ndips){
    vol.start<-volume(d) # get enclosure volume from depth
    v.sample<-dip.vol.flush(d) # how much did I scoop out?
    collected<-collected+dip.vol.flush(d) # track total collected so far
    d<-numeric.depth(new.depth(d.start=d,v.sample=v.sample)) # update depth
    # print(d)
  }
  return(collected)
}

sample.vol.hi<-function(depth,ndips=3){
  d<-depth
  # print(d)
  collected=0
  for(i in 1:ndips){
    vol.start<-volume(d) # get enclosure volume from depth
    v.sample<-dip.vol.hi(d) # how much did I scoop out?
    collected<-collected+dip.vol.hi(d) # track total collected so far
    d<-numeric.depth(new.depth(d.start=d,v.sample=v.sample)) # update depth
    # print(d)
  }
  return(collected)
}
