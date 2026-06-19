FROM mambaorg/micromamba:1.5.10

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PATH=/opt/conda/envs/emmaenv/bin:$PATH \
    DASH_HOST=0.0.0.0 \
    DASH_PORT=8050 \
    DASH_DEBUG=false

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

EXPOSE 8050

# Default: run the Dash dashboard on an externally reachable Docker interface.
# CLI pipeline override:
#   docker run --rm -it -p 8050:8050 -v "$PWD/input:/app/input" -v "$PWD/output:/app/output" emma-tool micromamba run -n emmaenv bash -c "cd code && python main.py"
CMD ["python", "dashboard/app.py"]
