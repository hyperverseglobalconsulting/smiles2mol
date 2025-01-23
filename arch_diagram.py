from diagrams import Diagram, Cluster
from diagrams.aws.network import APIGateway
from diagrams.aws.compute import Lambda
from diagrams.aws.storage import S3
from diagrams.onprem.ci import GithubActions
from diagrams.aws.network import Route53
from diagrams import Edge

# Create the architecture diagram
with Diagram("GitHub Pages to Lambda Architecture", show=False, outformat="jpg"):
    # Define the components
    github_pages = GithubActions("GitHub Pages")

    with Cluster("AWS Environment"):
        custom_domain = Route53("Custom Domain")
        api_gateway = APIGateway("API Gateway")
        lambda_function = Lambda("Lambda Function")
        s3_image_storage = S3("Image Storage")

    # Establish relationships
    github_pages >> api_gateway
    github_pages - Edge(style="dotted") >> custom_domain
    custom_domain - Edge(style="dotted") >> api_gateway
    api_gateway >> lambda_function
    lambda_function >> s3_image_storage
