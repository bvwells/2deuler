FROM gcc:latest

RUN apt-get update && apt-get install -y --no-install-recommends \ 
     cmake \ 
   && apt-get clean \ 
   && rm -rf /var/lib/apt/lists/* 

WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app
