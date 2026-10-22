######################## Animation function ####################################

#The function creates a gif animation from screenshots of each time steps 
#that are exported during the dynamic model run.

animate <- function(fps = 2){

## list file names and read in
imgs <- list.files("Outputs/animation/", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = fps)

## save to disk
image_write(image = img_animated,
            path = "Outputs/animation/model_animation.gif")

}