FROM us.gcr.io/broad-dsde-methods/glimpse:odelaneau_bd93ade

RUN apt-get update && apt-get install -y curl apt-transport-https ca-certificates gnupg && \
  echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
  apt-get update -y && apt-get install google-cloud-cli -y
