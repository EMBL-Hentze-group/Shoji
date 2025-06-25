
FROM python:3.12
COPY ./dist /tmp/
RUN cd /tmp/ \
    && ls *.whl| xargs -I whl pip install --no-cache-dir whl \
    && rm -rf *.whl \
    && pip cache purge
