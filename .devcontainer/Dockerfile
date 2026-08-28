
FROM mathworks/matlab:R2025b

USER matlab
WORKDIR /home/matlab

RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && HOME=${HOME} sudo ./mpm install \
    --release=R2025b --destination=/opt/matlab/R2025b \
    --products=Control_System_Toolbox Robust_Control_Toolbox Symbolic_Math_Toolbox

RUN bash -c "sudo apt-get update && sudo apt-get install --no-install-recommends -y git python3 python3-pip xvfb fluxbox octave octave-dev"

RUN bash -c "pipx upgrade matlab-proxy \
    && pipx inject --include-apps --include-deps matlab-proxy jupyter-matlab-proxy ipykernel jupyterlab octave_kernel"

USER root
ENV HOME=/home/matlab